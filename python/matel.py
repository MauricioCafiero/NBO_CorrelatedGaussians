"""
Matrix-element computation for NECG.

Analogous to matel.f90.  Computes the overlap (ov), kinetic-energy (tk),
electron-repulsion (ert), and electric-field (zht) contributions for a
pair of Correlated-Gaussian basis functions ii and jj under symmetry
term stf.

All indices are **1-based** on input (matching Fortran call sites) and
converted internally.

Key NumPy replacements:
  - Fortran ``matmul``        → ``@`` operator
  - Fortran DPOTRF + DPOTRI   → ``numpy.linalg.inv``
  - Kronecker products        → ``numpy.kron(A, np.eye(3))``
"""

from __future__ import annotations

import math

import numpy as np

import shared as sh
from utils import f0


def matel(stf: int, ii: int, jj: int) -> tuple[float, float, float, float]:
    """
    Compute matrix elements between basis functions *ii* and *jj* under
    symmetry term *stf* (all 1-based, matching Fortran convention).

    Returns
    -------
    ov  : overlap integral
    tk  : kinetic-energy integral
    ert : electron-repulsion integral (summed over all pairs)
    zht : electric-field (Stark) integral
    """
    # Convert to 0-based indices
    stf0 = stf - 1
    ii0 = ii - 1
    jj0 = jj - 1

    ne = sh.ne
    svh = sh.svh
    skp = sh.skp
    blk = svh + skp   # size of one basis-function block in iv

    # ------------------------------------------------------------------
    # Build lower-triangular Cholesky factors lk (bra) and ell (ket)
    # from the iv vector.
    # iv[(ii0*blk) .. (ii0*blk + svh - 1)] holds vech(lk)
    # iv[(jj0*blk) .. (jj0*blk + svh - 1)] holds vech(ell)
    # ------------------------------------------------------------------
    lk = np.zeros((ne, ne))
    ell = np.zeros((ne, ne))

    kk = ii0 * blk
    ll = jj0 * blk
    for i in range(ne):          # column (Fortran i)
        for j in range(i, ne):   # row   (Fortran j >= i → lower triangle)
            lk[j, i] = sh.iv[kk]
            ell[j, i] = sh.iv[ll]
            kk += 1
            ll += 1

    # Positive-semidefinite matrices  Ak = lk lk^T,  Al = ell ell^T
    ak = lk @ lk.T
    al = ell @ ell.T

    # ------------------------------------------------------------------
    # Apply symmetry operation on the ket:  al → sml^T · al · sml
    # Fortran: alt = al @ sml;  al = sml^T @ alt
    # ------------------------------------------------------------------
    sml_k = sh.sml[:, :, stf0]
    al = sml_k.T @ (al @ sml_k)

    # ------------------------------------------------------------------
    # Shift vectors  st1 (bra) and st2 (ket)
    # iv[(ii0*blk + svh) .. (ii0*blk + svh + skp - 1)]
    # ------------------------------------------------------------------
    st1 = sh.iv[ii0 * blk + svh: ii0 * blk + svh + skp].copy()
    st2 = sh.iv[jj0 * blk + svh: jj0 * blk + svh + skp].copy()

    # Apply symmetry on ket shift vector
    st2 = sh.sms[:, :, stf0] @ st2

    # ------------------------------------------------------------------
    # Matrix inverses via numpy (replacing DPOTRF + DPOTRI)
    # ------------------------------------------------------------------
    akl = ak + al
    akli = np.linalg.inv(akl)   # (ak + al)^{-1}

    # ------------------------------------------------------------------
    # Determinant factor.
    # Fortran: DPOTRF fills L (Cholesky of akl); dett = prod(L_ii)^2 = det(akl)
    # ------------------------------------------------------------------
    dett = abs(np.linalg.det(akl))   # = det(ak + al)

    # Product of diagonal elements of the Cholesky factors of ak and al
    # (det of lower-triangular lk / ell = product of diagonal)
    detak = math.sqrt(abs(np.prod(np.diag(lk)) ** 2))
    detal = math.sqrt(abs(np.prod(np.diag(ell)) ** 2))

    # ------------------------------------------------------------------
    # Kronecker products  (ne × ne → skp × skp block-diagonal matrices)
    # np.kron(A, eye(3)) produces the same result as the triple-k loop
    # in the Fortran source.
    # ------------------------------------------------------------------
    I3 = np.eye(3)
    akb = np.kron(ak, I3)
    alb = np.kron(al, I3)
    aklib = np.kron(akli, I3)
    massb = np.kron(sh.mass, I3)

    # ------------------------------------------------------------------
    # Overlap integral
    # i1 = st1^T akb st1,  i2 = st2^T alb st2
    # ee = akb·st1 + alb·st2
    # st = aklib · ee
    # i3 = ee^T · aklib · ee = ee^T · st
    # ov = exp(-(i1+i2)+i3) * sqrt(8^ne) * sqrt((detak·detal/dett)^3)
    # ------------------------------------------------------------------
    iv1 = akb @ st1
    i1 = float(st1 @ iv1)

    iv2 = alb @ st2
    i2 = float(st2 @ iv2)

    ee = iv1 + iv2
    st = aklib @ ee
    i3 = float(ee @ st)

    prefac = math.sqrt(8.0 ** ne) * math.sqrt((detak * detal / dett) ** 3)
    ov = math.exp(-(i1 + i2) + i3) * prefac

    # ------------------------------------------------------------------
    # Kinetic-energy integral
    # tk = 2·ov·(2·i3_ke + trace_ke)
    # i3_ke = (st - st1)^T · (akb · massb · alb) · (st - st2)
    # trace_ke = tr(massb · akb · aklib · alb)
    # ------------------------------------------------------------------
    im1 = massb @ alb
    im2 = akb @ im1
    dif1 = st - st1
    dif2 = st - st2
    ke_i3 = float(dif1 @ (im2 @ dif2))

    im1 = aklib @ alb
    im2 = akb @ im1
    im1 = massb @ im2
    ke_i2 = float(np.trace(im1))

    tk = 2.0 * ov * (2.0 * ke_i3 + ke_i2)

    # ------------------------------------------------------------------
    # Electron-repulsion and nuclear-attraction integrals
    # Loop over all pairs (pp, ll_) with 0 ≤ pp < ll_ ≤ ne.
    # pp = 0 → nucleus–electron;  pp ≥ 1 → electron–electron.
    # ------------------------------------------------------------------
    ert = 0.0
    # pre-compute ism1 and store for each pair
    for pp in range(ne + 1):
        for ll_ in range(pp + 1, ne + 1):
            # Build the J matrix (ne × ne) — Fortran 1-based indices
            jij = np.zeros((ne, ne))
            if pp == 0:
                # nucleus–electron: only diagonal element for particle ll_
                jij[ll_ - 1, ll_ - 1] = 1.0
            else:
                jij[pp - 1, pp - 1] = 1.0
                jij[ll_ - 1, ll_ - 1] = 1.0
                jij[pp - 1, ll_ - 1] = -1.0
                jij[ll_ - 1, pp - 1] = -1.0

            # ism1 = jij · akli,  jij ← akli · ism1 = akli · jij_orig · akli
            ism1 = jij @ akli
            jij_new = akli @ ism1

            # Kronecker product of updated jij
            jijb = np.kron(jij_new, I3)

            # i3_er = ee^T · jijb · ee
            iv3 = jijb @ ee
            i3_er = float(ee @ iv3)

            # i2_er = tr(ism1) = tr(jij_orig · akli)
            i2_er = float(np.trace(ism1))

            if i2_er == 0.0:
                continue

            i1_er = i3_er / i2_er
            er = (sh.ncv[pp] * sh.ncv[ll_]
                  * 2.0 * ov * f0(i1_er) / math.sqrt(i2_er * sh.PI))
            ert += er

    # ------------------------------------------------------------------
    # Electric-field (Stark) term
    # zht = ov · efs · sum_i( ncv[i] · st[3*(i-1)+2] )
    # where i runs over electrons (1-indexed), st[3*(i-1)+2] is the z-
    # component (index 3*i in Fortran 1-based = 3*(i-1)+2 in Python).
    # ------------------------------------------------------------------
    zht_sum = 0.0
    for i in range(1, ne + 1):          # Fortran i = 1..ne
        zht_sum += sh.ncv[i] * st[3 * i - 1]   # st(3*i) in Fortran 1-based
    zht = ov * zht_sum * sh.efs

    return ov, tk, ert, zht
