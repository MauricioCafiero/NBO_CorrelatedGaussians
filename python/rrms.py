"""
Dipole-moment matrix elements for NECG.

Analogous to rrms.f90.  Computes the overlap (ov) and a single Cartesian
component of the position expectation value (the n1-th element of the
aklib·ee vector) for basis-function pair (ii, jj) under symmetry terms
(stfe, stf).

All indices are **1-based** on input.
"""

from __future__ import annotations

import math

import numpy as np

import shared as sh


def rrms(stfe: int, stf: int, ii: int, jj: int, n1: int) -> tuple[float, float]:
    """
    Compute overlap and position-vector element for dipole calculation.

    Parameters
    ----------
    stfe : symmetry index for the bra (1-based)
    stf  : symmetry index for the ket (1-based)
    ii   : bra basis-function index (1-based)
    jj   : ket basis-function index (1-based)
    n1   : Cartesian component index (1-based, into skp-vector)

    Returns
    -------
    ov : overlap
    tk : n1-th component of aklib·ee, times ov
    """
    stfe0 = stfe - 1
    stf0 = stf - 1
    ii0 = ii - 1
    jj0 = jj - 1

    ne = sh.ne
    svh = sh.svh
    skp = sh.skp
    blk = svh + skp

    # ------------------------------------------------------------------
    # Build Cholesky factors
    # ------------------------------------------------------------------
    lk = np.zeros((ne, ne))
    ell = np.zeros((ne, ne))

    kk = ii0 * blk
    ll = jj0 * blk
    for i in range(ne):
        for j in range(i, ne):
            lk[j, i] = sh.iv[kk]
            ell[j, i] = sh.iv[ll]
            kk += 1
            ll += 1

    ak = lk @ lk.T
    al = ell @ ell.T

    # Apply symmetry on ket:  al → sml[:,:,stf]^T · al · sml[:,:,stf]
    sml_k = sh.sml[:, :, stf0]
    al = sml_k.T @ (al @ sml_k)

    # Apply symmetry on bra:  ak → sml[:,:,stfe]^T · ak · sml[:,:,stfe]
    sml_ke = sh.sml[:, :, stfe0]
    ak = sml_ke.T @ (ak @ sml_ke)

    # Shift vectors
    st1 = sh.iv[ii0 * blk + svh: ii0 * blk + svh + skp].copy()
    st2 = sh.iv[jj0 * blk + svh: jj0 * blk + svh + skp].copy()

    # Apply symmetry on ket shift
    st2 = sh.sms[:, :, stf0] @ st2
    # Apply symmetry on bra shift
    st1 = sh.sms[:, :, stfe0] @ st1

    # Matrix inverse
    akl = ak + al
    akli = np.linalg.inv(akl)

    # Determinants
    # Fortran: DPOTRF fills L (Cholesky of akl); dett = prod(L_ii)^2 = det(akl)
    dett = abs(np.linalg.det(akl))   # = det(ak + al)
    detak = math.sqrt(abs(np.prod(np.diag(lk)) ** 2))
    detal = math.sqrt(abs(np.prod(np.diag(ell)) ** 2))

    # Kronecker products
    I3 = np.eye(3)
    akb = np.kron(ak, I3)
    alb = np.kron(al, I3)
    aklib = np.kron(akli, I3)

    # ------------------------------------------------------------------
    # Overlap integral
    # ------------------------------------------------------------------
    iv1 = akb @ st1
    i1 = float(st1 @ iv1)

    iv2 = alb @ st2
    i2 = float(st2 @ iv2)

    ee = iv1 + iv2
    iv3 = aklib @ ee
    i3 = float(ee @ iv3)

    st = iv3   # st = aklib · ee

    prefac = math.sqrt(8.0 ** ne) * math.sqrt((detak * detal / dett) ** 3)
    ov = math.exp(-(i1 + i2) + i3) * prefac

    # n1-th component (1-based in Fortran, 0-based here)
    tk = st[n1 - 1] * ov

    return ov, tk
