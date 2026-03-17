"""
Utility functions for NECG calculations.

Provides:
  - biout(string)          : write to stdout and output file
  - get_current_time()     : return (date_str, time_str, zone_str)
  - vech(Ma, dim)          : extract lower-triangular elements of Ma
  - tab(ax, n, m)          : print matrix in tabular format
  - f0(z)                  : Boys function F0 via Pade approximation
  - df0(z)                 : derivative of F0
"""

from __future__ import annotations

import math
from datetime import datetime, timezone

import numpy as np

import shared as sh

# Path for the output file (mirrors unit-7 in the Fortran code)
_OUTPUT_FILE = "necg.out"


def biout(string: str) -> None:
    """Write *string* to stdout and append it to the output file."""
    print(string)
    with open(_OUTPUT_FILE, "a", encoding="utf-8") as fh:
        fh.write(string + "\n")


def get_current_time() -> tuple[str, str, str]:
    """Return (date, time, timezone) strings matching Fortran DATE_AND_TIME."""
    now = datetime.now(tz=timezone.utc)
    td = now.strftime("%Y%m%d")
    tt = now.strftime("%H%M%S") + ".000"
    tz = "+0000"
    return td, tt, tz


def vech(Ma: np.ndarray, dim: int) -> np.ndarray:
    """
    Return the half-vectorisation (vech) of the lower triangle of *Ma*.

    Elements are stored column by column (j = 0 .. dim-1), row from j
    downwards — matching the Fortran subroutine vech.f90.
    """
    vec = np.empty(dim * (dim + 1) // 2)
    k = 0
    for j in range(dim):
        for i in range(j, dim):
            vec[k] = Ma[i, j]
            k += 1
    return vec


def tab(ax: np.ndarray, n: int, m: int) -> None:
    """Print matrix *ax* (n rows × m cols) in tabular format (mirrors TAB in matrix_print.f90)."""
    num_blocks = m // 10
    for blk in range(num_blocks):
        jp = blk * 10
        jk = jp + 10
        header = "      " + "".join(f"{j + 1:24d}" for j in range(jp, jk))
        print(header)
        for i in range(n):
            row = f"{i + 1:3d}   " + "".join(f"{ax[i, j]:23.15f}" for j in range(jp, jk))
            print(row)
    ma = num_blocks * 10
    if ma < m:
        header = "      " + "".join(f"{j + 1:24d}" for j in range(ma, m))
        print(header)
        for i in range(n):
            row = f"{i + 1:3d}   " + "".join(f"{ax[i, j]:23.15f}" for j in range(ma, m))
            print(row)


def f0(z: float) -> float:
    """
    Boys function F0(z) via Pade approximation (error_function_approximation.f90).

    Uses Pade coefficients stored in *shared*.  Falls back to the
    asymptotic form sqrt(pi/(2z)) / sqrt(2) for large z.
    """
    if z <= 16.3578:
        num = (((((sh.C1 * z + sh.C2) * z + sh.C3) * z + sh.C4) * z + sh.C5) * z + 1.0)
        den = ((((((sh.C6 * z + sh.C7) * z + sh.C8) * z + sh.C9) * z + sh.C10) * z + sh.C11) * z + 1.0)
        return math.sqrt(num / den)
    r = 2.0 * z
    return math.sqrt(math.pi / (2.0 * r))


def df0(z: float) -> float:
    """Derivative of Boys function F0 (derror.f90)."""
    if z <= 16.3578:
        num = (((((sh.C1 * z + sh.C2) * z + sh.C3) * z + sh.C4) * z + sh.C5) * z + 1.0)
        den = ((((((sh.C6 * z + sh.C7) * z + sh.C8) * z + sh.C9) * z + sh.C10) * z + sh.C11) * z + 1.0)
        f0_val = math.sqrt(num / den)
        anum = sh.C1 * z ** 5 + sh.C2 * z ** 4 + sh.C3 * z ** 3 + sh.C4 * z ** 2 + sh.C5 * z + 1.0
        aden = sh.C6 * z ** 6 + sh.C7 * z ** 5 + sh.C8 * z ** 4 + sh.C9 * z ** 3 + sh.C10 * z ** 2 + sh.C11 * z + 1.0
        danum = 5.0 * sh.C1 * z ** 4 + 4.0 * sh.C2 * z ** 3 + 3.0 * sh.C3 * z ** 2 + 2.0 * sh.C4 * z + sh.C5
        daden = 6.0 * sh.C6 * z ** 5 + 5.0 * sh.C7 * z ** 4 + 4.0 * sh.C8 * z ** 3 + 3.0 * sh.C9 * z ** 2 + 2.0 * sh.C10 * z + sh.C11
        return (1.0 / (2.0 * f0_val)) * ((aden * danum - anum * daden) / aden ** 2)
    # Asymptotic branch: d/dz sqrt(pi/(4z)) = -sqrt(pi) / (4 * z^(3/2)) / 2 ...
    # Simplifies to: -pi / (8 * z^2) / f0(z)
    return -math.sqrt(sh.PI / (4.0 * z)) / (2.0 * z)
