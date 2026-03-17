"""
Objective function and gradient for the NECG optimiser.

Analogous to outprod.f90.  Given a parameter vector *x* (concatenation of
iv and C), recomputes H and S, evaluates the Rayleigh quotient f = <C|H|C>
/ <C|S|C>, and builds the analytical gradient g.

MPI calls removed; runs on a single process.
"""

from __future__ import annotations

import numpy as np

import shared as sh
from matel import matel
from make_grad import make_grad


def outprod(x: np.ndarray) -> tuple[float, np.ndarray]:
    """
    Evaluate energy and gradient for parameter vector *x*.

    Parameters
    ----------
    x : 1-D array of length siv + M
        First siv elements → nonlinear parameters iv.
        Next M elements    → linear coefficients C.

    Returns
    -------
    f : float
        Rayleigh quotient (energy).
    g : ndarray, shape (siv+M,)
        Analytical gradient of f w.r.t. x.
    """
    # ------------------------------------------------------------------
    # Unpack x into iv and C
    # ------------------------------------------------------------------
    sh.iv[:] = x[:sh.siv]
    sh.C[:] = x[sh.siv:]

    # ------------------------------------------------------------------
    # Build H and S matrices
    # ------------------------------------------------------------------
    sh.H[:] = 0.0
    sh.S[:] = 0.0

    for i in range(1, sh.M + 1):
        for j in range(1, i + 1):
            for stf in range(1, sh.nst + 1):
                ov, tk, ert, zht = matel(stf, i, j)
                sh.H[i - 1, j - 1] += sh.scvv[stf - 1] * (tk + ert - zht)
                sh.S[i - 1, j - 1] += sh.scvv[stf - 1] * ov

    # Fill upper triangle by symmetry
    for j in range(sh.M):
        for i in range(j):
            sh.H[i, j] = sh.H[j, i]
            sh.S[i, j] = sh.S[j, i]

    # ------------------------------------------------------------------
    # Rayleigh quotient  f = C^T H C / C^T S C
    # ------------------------------------------------------------------
    cHc = float(sh.C @ sh.H @ sh.C)
    cSc = float(sh.C @ sh.S @ sh.C)
    f = cHc / cSc
    norm = 1.0 / cSc

    # ------------------------------------------------------------------
    # Analytical gradient
    # ------------------------------------------------------------------
    make_grad(f, norm)

    g = sh.grad.copy()
    return f, g
