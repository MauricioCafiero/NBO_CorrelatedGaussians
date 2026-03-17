"""
N-electron Correlated Gaussians (NECG) — Python implementation.

Project title : N-particle Non-Born-Oppenheimer Correlated Gaussians
Original code : ADAMOWICZ RESEARCH [UofA chemistry]
Original prog : Mauricio Cafiero
Python port   : See repository README

Analogue of necg.f90.  MPI parallelism has been removed; numpy / scipy
are used for all linear-algebra operations.

Usage
-----
    python necg.py [input_file]

The input file defaults to ``oldquit.dat`` (matching the Fortran version).
Output is written to stdout **and** to ``necg.out``.
"""

from __future__ import annotations

import sys

import numpy as np
from scipy import linalg
from scipy.optimize import minimize

import shared as sh
from read_input import read_input
from matel import matel
from outprod import outprod
from rrms import rrms
from utils import biout, get_current_time


def _build_matrices() -> None:
    """Fill H, S, Kinetic and Potential from scratch using current iv / C."""
    sh.H[:] = 0.0
    sh.S[:] = 0.0
    sh.Kinetic[:] = 0.0
    sh.Potential[:] = 0.0

    for i in range(1, sh.M + 1):
        for j in range(i, sh.M + 1):
            for stf in range(1, sh.nst + 1):
                ov, tk, ert, zht = matel(stf, j, i)
                sh.Kinetic[j - 1, i - 1] += sh.scvv[stf - 1] * tk
                sh.Potential[j - 1, i - 1] += sh.scvv[stf - 1] * ert
                sh.H[j - 1, i - 1] += sh.scvv[stf - 1] * (tk + ert - zht)
                sh.S[j - 1, i - 1] += sh.scvv[stf - 1] * ov

    # Symmetrise
    for j in range(sh.M):
        for i in range(j):
            sh.H[i, j] = sh.H[j, i]
            sh.S[i, j] = sh.S[j, i]
            sh.Kinetic[i, j] = sh.Kinetic[j, i]
            sh.Potential[i, j] = sh.Potential[j, i]


def _rayleigh_quotient() -> tuple[float, float, float, float]:
    """Return (energy, kinetic, potential, virial) using current C, H, S."""
    cHc = float(sh.C @ sh.H @ sh.C)
    cSc = float(sh.C @ sh.S @ sh.C)
    cKc = float(sh.C @ sh.Kinetic @ sh.C)
    cPc = float(sh.C @ sh.Potential @ sh.C)
    energy = cHc / cSc
    kin = cKc / cSc
    pot = cPc / cSc
    virial = -2.0 * kin / pot
    return energy, kin, pot, virial


def _print_energetics(energy: float, kin: float, pot: float, virial: float) -> None:
    print("*" * 69)
    print()
    print(f"Rayleigh Quotient: {energy:25.15f}")
    print()
    print(f"Average Kinetic  : {kin:25.15f}")
    print()
    print(f"Average Potential: {pot:25.15f}")
    print()
    print(f"Virial Coef.     : {virial:25.15f}")
    print()
    print("*" * 69)


def main(input_file: str = "oldquit.dat") -> None:
    """Main driver — mirrors the flow of necg.f90."""

    # ------------------------------------------------------------------
    # Initialise output
    # ------------------------------------------------------------------
    import os
    if os.path.exists("necg.out"):
        os.remove("necg.out")

    biout("N-particle Non-Born-Oppenheimer Correlated Gaussians")
    biout("=" * 68)
    biout("ADAMOWICZ RESEARCH                  [UofA chemistry]")
    biout("Programming: Mauricio Cafiero")
    biout("=" * 68)
    biout(":PROGRAM BEGINS:")
    td, tt, tz = get_current_time()
    biout(f"Date: {td} | Time: {tt} / {tz} gmt")

    # ------------------------------------------------------------------
    # Read input
    # ------------------------------------------------------------------
    read_input(input_file)

    print(f"M, NE, NST = {sh.M}, {sh.ne}, {sh.nst}")
    print(f"Electric Field = {sh.efs}")
    print("OUTPUT FILENAME = necg.out")

    # ------------------------------------------------------------------
    # Build initial H and S matrices
    # ------------------------------------------------------------------
    import time
    t0 = time.time()
    _build_matrices()
    t1 = time.time()
    print(f"Time for matrix calc. is {t1 - t0:.3f} seconds")

    # ------------------------------------------------------------------
    # Solve generalised eigenvalue problem  H·v = λ·S·v
    # to initialise C (scipy.linalg.eigh replaces DSYGV)
    # ------------------------------------------------------------------
    try:
        eigenvalues, eigenvectors = linalg.eigh(sh.H, sh.S)
        print(f"Eigenvalue = {eigenvalues[0]:.15f}")
        sh.C[:] = eigenvectors[:, 0]
    except linalg.LinAlgError as exc:
        print(f"Warning: generalised eigenvalue problem failed: {exc}")

    print()
    print("***** The Input Vectors *****")
    print(sh.C)
    print("***** End Input Vectors *****")
    print()

    # ------------------------------------------------------------------
    # Initial energetics
    # ------------------------------------------------------------------
    energy, kin, pot, virial = _rayleigh_quotient()
    _print_energetics(energy, kin, pot, virial)

    # ------------------------------------------------------------------
    # Assemble initial parameter vector for optimiser
    # ------------------------------------------------------------------
    hiv = np.concatenate([sh.iv.copy(), sh.C.copy()])

    # ------------------------------------------------------------------
    # Optimisation using L-BFGS-B (scipy analogue of LMQN)
    # scipy.optimize.minimize(method="L-BFGS-B") provides a limited-
    # memory quasi-Newton method analogous to the LMQN Fortran routine.
    # ------------------------------------------------------------------
    t0 = time.time()
    result = minimize(
        outprod,
        hiv,
        method="L-BFGS-B",
        jac=True,                           # outprod returns (f, g)
        options={
            "maxiter": 20,
            "maxfun": 150 * (sh.siv + sh.M),
            "ftol": 2.220446649250313e-16,
            "gtol": np.sqrt(2.220446649250313e-16),
            "iprint": 1,
        },
    )
    t1 = time.time()
    print(f"Time for energy/gradient calc. is {t1 - t0:.3f} seconds")
    print(f"Energy is {result.fun:.15f} hartree")
    print(f"Optimiser message: {result.message}")

    # ------------------------------------------------------------------
    # Unpack optimised parameters
    # ------------------------------------------------------------------
    sh.iv[:] = result.x[:sh.siv]
    sh.C[:] = result.x[sh.siv:]

    print()
    print("***** The Output Vectors *****")
    print(sh.C)
    print()
    print("***** End Output Vectors *****")
    print()

    # ------------------------------------------------------------------
    # Final matrix build with optimised parameters
    # ------------------------------------------------------------------
    _build_matrices()
    energy, kin, pot, virial = _rayleigh_quotient()
    _print_energetics(energy, kin, pot, virial)

    # ------------------------------------------------------------------
    # Dipole moment calculation (replaces the rrms loop in necg.f90)
    # ------------------------------------------------------------------
    d7 = d8 = d9 = 0.0
    kk = 1
    for ii in range(1, sh.ne + 1):
        for jj_comp in range(1, 4):
            h_dip = np.zeros((sh.M, sh.M))
            s_dip = np.zeros((sh.M, sh.M))

            for stfe in range(1, sh.nst + 1):
                for stf in range(1, sh.nst + 1):
                    for i in range(1, sh.M + 1):
                        for j in range(i, sh.M + 1):
                            ov_r, tk_r = rrms(stfe, stf, j, i, kk)
                            h_dip[j - 1, i - 1] += sh.scvv[stf - 1] * sh.scvv[stfe - 1] * tk_r
                            s_dip[j - 1, i - 1] += sh.scvv[stf - 1] * sh.scvv[stfe - 1] * ov_r

            # Symmetrise
            for j in range(sh.M):
                for i in range(j):
                    h_dip[i, j] = h_dip[j, i]
                    s_dip[i, j] = s_dip[j, i]

            cHc_d = float(sh.C @ h_dip @ sh.C)
            cSc_d = float(sh.C @ s_dip @ sh.C)
            d3 = cHc_d / cSc_d if cSc_d != 0.0 else 0.0

            if jj_comp == 1:
                d7 += sh.ncv[ii] * d3
            elif jj_comp == 2:
                d8 += sh.ncv[ii] * d3
            else:
                d9 += sh.ncv[ii] * d3

            kk += 1

    DEBYE = 2.541765e0
    d4 = d7 * DEBYE
    d5 = d8 * DEBYE
    d6 = d9 * DEBYE

    print()
    print(f"Dipole x-component = {d7:.15f}  ( {d4:.6f} Debye )")
    print()
    print(f"Dipole y-component = {d8:.15f}  ( {d5:.6f} Debye )")
    print()
    print(f"Dipole z-component = {d9:.15f}  ( {d6:.6f} Debye )")
    print()
    print("*" * 69)

    biout("PROGRAM TERMINATES")
    td, tt, tz = get_current_time()
    biout(f"Date: {td} | Time: {tt} / {tz} gmt")


if __name__ == "__main__":
    input_file = sys.argv[1] if len(sys.argv) > 1 else "oldquit.dat"
    main(input_file)
