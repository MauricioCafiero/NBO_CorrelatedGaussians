"""
Gradient assembly for NECG.

Analogous to make_grad.f90.  Builds the full gradient vector ``sh.grad``
from the individual element gradients computed by gradel().

MPI calls removed; the loops run sequentially on a single process.
"""

from __future__ import annotations

import numpy as np

import shared as sh
from gradel import gradel


def make_grad(eng: float, norm: float) -> None:
    """
    Build the full gradient vector and store it in ``shared.grad``.

    Parameters
    ----------
    eng  : current energy (Rayleigh quotient)
    norm : 1 / <C|S|C>   (normalisation factor)
    """
    svh = sh.svh
    skp = sh.skp
    blk = svh + skp     # size of one parameter block

    # ------------------------------------------------------------------
    # Accumulate element gradients into DH and DS
    # ------------------------------------------------------------------
    sh.DH = np.zeros((sh.siv, sh.M))
    sh.DS = np.zeros((sh.siv, sh.M))

    for l in range(1, sh.M + 1):
        j = 0
        for k in range(1, sh.M + 1):
            tgvt = np.zeros(blk)
            tgvst = np.zeros(blk)

            for stf in range(1, sh.nst + 1):
                tgv, tgvs = gradel(stf, l, k)
                tgvt += sh.scvv[stf - 1] * tgv
                tgvst += sh.scvv[stf - 1] * tgvs

            # Store into DH / DS
            # Fortran: dh((l-k)*(svh+skp)+j .. , k) = tgvt
            # In Python the column index for H(l,k) is k-1.
            row_start = (l - k) * blk + j
            for i in range(blk):
                sh.DH[row_start + i, k - 1] = tgvt[i]
                sh.DS[row_start + i, k - 1] = tgvst[i]
            j += blk

    # ------------------------------------------------------------------
    # Build GRADM: (M × siv)
    # GRADM[i, j] = norm * (DH[j, i] - eng * DS[j, i])
    # ------------------------------------------------------------------
    GRADM = norm * (sh.DH - eng * sh.DS).T   # shape (M, siv)

    # ------------------------------------------------------------------
    # Non-linear parameter gradient (GRADA): length siv
    # GRADA[ic] = sum_i c_i * (2·c_j·GRADM[i,ic]  if i≠j
    #                          or c_i·c_j·GRADM[i,ic]  if i=j)
    # ------------------------------------------------------------------
    ii_blk = blk
    GRADA = np.zeros(sh.siv)
    for j in range(1, sh.M + 1):
        for k in range(1, ii_blk + 1):
            ic = (j - 1) * ii_blk + k - 1
            val = 0.0
            for i in range(1, sh.M + 1):
                if i != j:
                    val += 2.0 * sh.C[i - 1] * sh.C[j - 1] * GRADM[i - 1, ic]
                else:
                    val += sh.C[i - 1] * sh.C[j - 1] * GRADM[i - 1, ic]
            GRADA[ic] = val

    # ------------------------------------------------------------------
    # Linear coefficient gradient (GRADC): length M
    # GRADC[i] = 2·norm · sum_j (H[i,j] - eng·S[i,j]) · C[j]
    # ------------------------------------------------------------------
    GRADC = 2.0 * norm * (sh.H - eng * sh.S) @ sh.C

    # ------------------------------------------------------------------
    # Assemble full gradient
    # ------------------------------------------------------------------
    sh.grad = np.concatenate([GRADA, GRADC])
