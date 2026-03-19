#!/usr/bin/env python3
"""
Generate an input file for the LiH molecule with any number of basis functions.

The reference file ``oldquit.dat`` (12-function LiH basis) is parsed to
extract physical parameters (mass tensor, charges, symmetry matrices) and the
statistical properties of the nonlinear parameter blocks.  A new file is then
written with a user-specified number of basis functions *M*, using novel but
statistically similar Cholesky and shift parameters sampled from distributions
that match the reference data.

Usage
-----
    python generate_lih_input.py --M 16 --output lih_16.dat
    python generate_lih_input.py --M 8  --output lih_8.dat  --seed 42
    python generate_lih_input.py --M 20  # writes to lih_M20.dat

The generated file can be used directly as input to the NECG code
(``python/necg.py``).
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

# Default path to the reference input file (one level above this script).
_HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_REFERENCE = os.path.join(_HERE, "..", "oldquit.dat")


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def _parse_reference(filename: str) -> dict:
    """Parse *filename* and return all parameters as a dictionary.

    The returned dictionary contains:

    ``M_ref`` : int
        Number of basis functions in the reference file.
    ``ne`` : int
        Number of pseudo-particles (electrons + nuclei, non-BO).
    ``nst`` : int
        Number of symmetry terms.
    ``efs`` : float
        Electric-field strength.
    ``svh`` : int
        ``ne*(ne+1)//2`` — number of lower-triangular Cholesky elements.
    ``skp`` : int
        ``3*ne`` — length of the shift vector per basis function.
    ``mass`` : ndarray, shape (ne, ne)
        Inverse-mass tensor.
    ``ncv`` : ndarray, shape (ne+1,)
        Charge values.
    ``iv_blocks`` : ndarray, shape (M_ref, svh+skp)
        Nonlinear parameters, one row per reference basis function.
    ``C_ref`` : ndarray, shape (M_ref,)
        Linear expansion coefficients from the reference.
    ``scvv`` : ndarray, shape (nst,)
        Symmetry-term weights (raw, before /nst normalisation used by NECG).
    ``sml`` : ndarray, shape (nst, ne, ne)
        Symmetry permutation matrices for the Gaussian exponent matrices.
    ``sms`` : ndarray, shape (nst, skp, skp)
        Symmetry permutation matrices for the shift vectors.
    """
    try:
        with open(filename, "r", encoding="utf-8") as fh:
            tokens = fh.read().split()
    except OSError as exc:
        sys.exit(f"ERROR: Cannot open reference file '{filename}': {exc}")

    pos = 0

    def _int() -> int:
        nonlocal pos
        v = int(tokens[pos]); pos += 1; return v

    def _flt() -> float:
        nonlocal pos
        v = float(tokens[pos].lower().replace("d", "e")); pos += 1; return v

    M_ref = _int()
    ne    = _int()
    nst   = _int()
    efs   = _flt()

    skp = 3 * ne
    svh = ne * (ne + 1) // 2
    siv = M_ref * (svh + skp)

    mass = np.array([[_flt() for _ in range(ne)] for _ in range(ne)])
    ncv  = np.array([_flt() for _ in range(ne + 1)])
    iv   = np.array([_flt() for _ in range(siv)])
    C_ref = np.array([_flt() for _ in range(M_ref)])
    scvv = np.array([_flt() for _ in range(nst)])

    sml = np.zeros((nst, ne, ne))
    for k in range(nst):
        for i in range(ne):
            for j in range(ne):
                sml[k, i, j] = _flt()

    sms = np.zeros((nst, skp, skp))
    for k in range(nst):
        for i in range(skp):
            for j in range(skp):
                sms[k, i, j] = _flt()

    return {
        "M_ref": M_ref,
        "ne": ne,
        "nst": nst,
        "efs": efs,
        "svh": svh,
        "skp": skp,
        "mass": mass,
        "ncv": ncv,
        "iv_blocks": iv.reshape(M_ref, svh + skp),
        "C_ref": C_ref,
        "scvv": scvv,
        "sml": sml,
        "sms": sms,
    }


# ---------------------------------------------------------------------------
# Basis-set generation
# ---------------------------------------------------------------------------

def _diagonal_positions(ne: int, svh: int) -> list[int]:
    """Return the indices of diagonal Cholesky elements in the packed vector.

    The lower-triangular Cholesky factor L (ne×ne) is stored column-by-column:
    L[0,0], L[1,0], ..., L[ne-1,0], L[1,1], L[2,1], ..., L[ne-1,1], ...

    Diagonal element L[i,i] occupies index  i*ne - i*(i-1)//2  in the
    packed vector.
    """
    diag = []
    idx = 0
    for col in range(ne):
        diag.append(idx)       # L[col, col] is the first element of column col
        idx += (ne - col)      # column col has (ne - col) elements
    return diag                # length == ne


def _fit_lognormal(mean: float, std: float) -> tuple[float, float]:
    """Return (mu, sigma) of the log-normal distribution that reproduces *mean*
    and *std* for the underlying normal variate.

    If mean <= 0 or std is very small, fall back to a safe default.
    """
    if mean <= 0:
        mean = abs(mean) + 1e-3
    if std < 1e-10:
        std = mean * 0.05
    cv2 = (std / mean) ** 2
    sigma2 = np.log1p(cv2)
    mu = np.log(mean) - 0.5 * sigma2
    return mu, np.sqrt(sigma2)


def _generate_iv_block(
    col_mean: np.ndarray,
    col_std:  np.ndarray,
    diag_idx: list[int],
    rng:      np.random.Generator,
) -> np.ndarray:
    """Sample one basis-function parameter block (svh + skp elements).

    * **Off-diagonal Cholesky elements** and **shift-vector components**:
      drawn from ``Normal(mean, std)`` for each column.
    * **Diagonal Cholesky elements** (must be strictly positive):
      drawn from a ``LogNormal`` distribution fitted to the reference mean
      and std, which guarantees positivity.
    """
    block = rng.normal(col_mean, col_std)

    # Replace diagonal Cholesky elements with log-normal samples.
    for idx in diag_idx:
        mu, sigma = _fit_lognormal(col_mean[idx], col_std[idx])
        block[idx] = rng.lognormal(mu, sigma)

    return block


# ---------------------------------------------------------------------------
# File formatting
# ---------------------------------------------------------------------------

def _fmt_float(v: float) -> str:
    """Format a float in the style of the original Fortran input files.

    Values very close to zero or exactly representable as short decimals are
    written compactly; others use 15-digit scientific notation.
    """
    if v == 0.0:
        return "0.0"
    return f"{v:25.15E}"


def _write_matrix_row(values: np.ndarray) -> str:
    """Return a single row of a matrix as a space-separated string."""
    return "  ".join(_fmt_float(v) for v in values)


def _write_lih_input(params: dict, M_new: int, rng: np.random.Generator) -> str:
    """Build and return the complete input file as a string.

    Parameters
    ----------
    params :
        Dictionary returned by :func:`_parse_reference`.
    M_new :
        Desired number of basis functions.
    rng :
        NumPy random Generator (seeded or default).
    """
    ne  = params["ne"]
    nst = params["nst"]
    efs = params["efs"]
    svh = params["svh"]
    skp = params["skp"]
    iv_blocks = params["iv_blocks"]   # (M_ref, svh+skp)

    # Per-column statistics from reference basis functions.
    col_mean = np.mean(iv_blocks, axis=0)
    col_std  = np.std(iv_blocks,  axis=0)
    # Prevent degenerate distributions where std ≈ 0 (e.g. near-zero columns).
    col_std  = np.where(col_std < 1e-10, 1e-6, col_std)

    diag_idx = _diagonal_positions(ne, svh)

    # Generate M_new basis-function parameter blocks.
    new_iv_blocks = np.stack(
        [_generate_iv_block(col_mean, col_std, diag_idx, rng) for _ in range(M_new)]
    )

    # Linear coefficients: equal-weight initialisation.
    C_new = np.full(M_new, 1.0 / M_new)

    # ------------------------------------------------------------------ #
    # Assemble output lines                                               #
    # ------------------------------------------------------------------ #
    lines: list[str] = []

    # --- Header scalars -------------------------------------------------
    lines.append(str(M_new))
    lines.append(str(ne))
    lines.append(str(nst))
    lines.append(f"{efs:.10E}".replace("E", "d"))
    lines.append("")

    # --- Inverse-mass tensor  (ne × ne) ---------------------------------
    for row in params["mass"]:
        lines.append(_write_matrix_row(row))
    lines.append("")

    # --- Charge values  (ne+1 values, one per line) ---------------------
    for v in params["ncv"]:
        lines.append(f"{v:.1f}")
    lines.append("")
    lines.append("")

    # --- Nonlinear parameters: M_new blocks, each of length svh+skp -----
    for k in range(M_new):
        # Cholesky elements (svh values)
        for v in new_iv_blocks[k, :svh]:
            lines.append(f"   {_fmt_float(v)}")
        # Shift-vector elements (skp values)
        for v in new_iv_blocks[k, svh:]:
            lines.append(f"   {_fmt_float(v)}")
        lines.append("")

    # --- Linear expansion coefficients  (M_new values) ------------------
    lines.append("")
    for v in C_new:
        lines.append(f"  {_fmt_float(v)}")
    lines.append("")

    # --- Symmetry weights  scvv (nst values) ----------------------------
    for v in params["scvv"]:
        # Write as integer if the value is whole, otherwise as float.
        if v == int(v):
            lines.append(str(int(v)))
        else:
            lines.append(f"{v:.6g}")
    lines.append("")

    # --- Symmetry matrices sml  (nst × ne × ne) -------------------------
    for k in range(nst):
        for i in range(ne):
            lines.append("  ".join(f"{v:.1f}" for v in params["sml"][k, i, :]))
        lines.append("")

    # --- Symmetry matrices sms  (nst × skp × skp) -----------------------
    for k in range(nst):
        for i in range(skp):
            lines.append("  ".join(f"{v:.1f}" for v in params["sms"][k, i, :]))
        lines.append("")

    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def generate_lih_input(
    M: int,
    reference: str = DEFAULT_REFERENCE,
    seed: int | None = None,
    output: str | None = None,
) -> str:
    """Generate an LiH input file with *M* basis functions.

    Parameters
    ----------
    M :
        Number of correlated Gaussian basis functions (must be >= 1).
    reference :
        Path to the reference ``oldquit.dat`` file.  Defaults to
        ``../oldquit.dat`` relative to this script.
    seed :
        Optional integer random seed for reproducibility.
    output :
        If given, the generated file is written to this path.  The
        content is also returned as a string regardless.

    Returns
    -------
    str
        The complete text of the generated input file.
    """
    if M < 1:
        raise ValueError(f"M must be >= 1, got {M}")

    params = _parse_reference(reference)
    rng    = np.random.default_rng(seed)
    text   = _write_lih_input(params, M, rng)

    if output:
        with open(output, "w", encoding="utf-8") as fh:
            fh.write(text)
        print(f"Generated LiH input: M={M} basis functions → '{output}'")

    return text


# ---------------------------------------------------------------------------
# Self-test (no external tools required)
# ---------------------------------------------------------------------------

def _self_test(reference: str = DEFAULT_REFERENCE) -> None:
    """Minimal round-trip validation.

    1. Generate a file for M=8 and M=16 with a fixed seed.
    2. Re-parse each generated file with :func:`_parse_reference`.
    3. Assert that the header scalars and structural dimensions are correct.
    4. Assert that all diagonal Cholesky elements are strictly positive.
    5. Assert that the symmetry matrices are reproduced exactly.
    """
    import tempfile

    ref = _parse_reference(reference)
    ne  = ref["ne"]
    nst = ref["nst"]
    svh = ref["svh"]
    skp = ref["skp"]
    diag_idx = _diagonal_positions(ne, svh)

    for M_test in (4, 8, 16):
        rng = np.random.default_rng(0)
        text = _write_lih_input(ref, M_test, rng)

        # Write to a temporary file and re-parse it.
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".dat", delete=False, encoding="utf-8"
        ) as tmp:
            tmp.write(text)
            tmp_path = tmp.name

        try:
            parsed = _parse_reference(tmp_path)
        finally:
            os.unlink(tmp_path)

        # --- Structural checks ------------------------------------------
        assert parsed["M_ref"] == M_test,  f"M mismatch: {parsed['M_ref']} != {M_test}"
        assert parsed["ne"]    == ne,      f"ne mismatch"
        assert parsed["nst"]   == nst,     f"nst mismatch"
        assert parsed["svh"]   == svh,     f"svh mismatch"
        assert parsed["skp"]   == skp,     f"skp mismatch"
        assert parsed["iv_blocks"].shape == (M_test, svh + skp)
        assert parsed["C_ref"].shape     == (M_test,)

        # --- Positivity of diagonal Cholesky elements -------------------
        chol = parsed["iv_blocks"][:, :svh]
        for idx in diag_idx:
            assert np.all(chol[:, idx] > 0), (
                f"Diagonal Cholesky element at position {idx} is non-positive"
            )

        # --- Symmetry matrices reproduced exactly -----------------------
        assert np.allclose(parsed["sml"], ref["sml"]), "sml mismatch"
        assert np.allclose(parsed["sms"], ref["sms"]), "sms mismatch"
        assert np.allclose(parsed["ncv"], ref["ncv"]), "ncv mismatch"

        print(f"  PASS  M={M_test:3d}: shape OK, diagonals positive, "
              f"symmetry matrices exact")

    print("All self-tests passed.")


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Generate an NECG input file for the LiH molecule with a "
            "user-specified number of correlated Gaussian basis functions."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--M", "-M",
        type=int,
        default=None,
        help="Number of basis functions to generate (required unless --test is used).",
    )
    p.add_argument(
        "--output", "-o",
        default=None,
        help=(
            "Output file path.  If omitted, the file is written to "
            "'lih_M<M>.dat' in the current working directory."
        ),
    )
    p.add_argument(
        "--reference", "-r",
        default=DEFAULT_REFERENCE,
        help="Path to the reference LiH input file (oldquit.dat).",
    )
    p.add_argument(
        "--seed", "-s",
        type=int,
        default=None,
        help="Integer random seed for reproducibility (optional).",
    )
    p.add_argument(
        "--test",
        action="store_true",
        help="Run the built-in self-tests and exit.",
    )
    return p


def main(argv: list[str] | None = None) -> None:
    parser = _build_parser()
    args   = parser.parse_args(argv)

    if args.test:
        print("Running self-tests …")
        _self_test(args.reference)
        return

    if args.M is None:
        parser.error("--M is required when not using --test")

    output = args.output or f"lih_M{args.M}.dat"
    generate_lih_input(
        M=args.M,
        reference=args.reference,
        seed=args.seed,
        output=output,
    )


if __name__ == "__main__":
    main()
