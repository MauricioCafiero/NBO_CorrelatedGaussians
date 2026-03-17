"""
Gradient of matrix elements for NECG.

Analogous to gradel.f90.  Computes analytical gradients of the overlap,
kinetic-energy, and electron-repulsion integrals with respect to the
nonlinear parameters (Cholesky factors and shift vectors).

All indices are **1-based** on input (matching Fortran call sites).

Returns
-------
tgvl  : gradient of the Hamiltonian element H(ii,jj)
tgvls : gradient of the overlap element S(ii,jj)
Both are 1-D arrays of length svh+skp.
"""

from __future__ import annotations

import math

import numpy as np

import shared as sh
from utils import f0, df0, vech


def gradel(stf: int, ii: int, jj: int) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute analytical gradients of matrix elements for basis-function
    pair (*ii*, *jj*) under symmetry term *stf* (all 1-based).

    Returns ``(tgvl, tgvls)`` — both are length-(svh+skp) arrays.
    """
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

    # Apply symmetry on ket:  al → sml^T · al · sml
    sml_k = sh.sml[:, :, stf0]
    al = sml_k.T @ (al @ sml_k)

    # Shift vectors
    st1 = sh.iv[ii0 * blk + svh: ii0 * blk + svh + skp].copy()
    st2 = sh.iv[jj0 * blk + svh: jj0 * blk + svh + skp].copy()
    st2 = sh.sms[:, :, stf0] @ st2

    # Matrix inverses
    akl = ak + al
    akli = np.linalg.inv(akl)
    aki = np.linalg.inv(ak)

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
    massb = np.kron(sh.mass, I3)

    # ------------------------------------------------------------------
    # Overlap integral
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
    # Gradient of S w.r.t. shift sk (bra shift, 'ogs')
    # ogs = 2·ov·akb·(st - st1)
    # ------------------------------------------------------------------
    ogs = 2.0 * ov * (akb @ (st - st1))

    # ------------------------------------------------------------------
    # Gradient of S w.r.t. vech(lk)  — 'ogl'
    # Several terms built via vech()
    # ------------------------------------------------------------------
    ogl = np.zeros(svh)

    ism1 = aki @ lk
    ogl += 1.5 * vech(ism1, ne)

    ism1 = akli @ lk
    ogl -= 3.0 * vech(ism1, ne)

    # Reshape st1 into (3, ne) column form — matching Fortran dot1 layout
    # Fortran stores as dot1(j,i) = st1(3*(i-1)+j) for j=1..3, i=1..ne
    dot1 = st1.reshape(ne, 3).T   # shape (3, ne)
    dot2 = dot1.T                 # shape (ne, 3)
    ism1 = dot2 @ dot1            # (ne, ne)
    ism2 = ism1 @ lk
    ogl -= 2.0 * vech(ism2, ne)

    dot1_ee = ee.reshape(ne, 3).T  # shape (3, ne)
    dot2_ee = dot1_ee.T            # shape (ne, 3)
    ism1 = dot2_ee @ dot1_ee
    ism2 = ism1 @ akli
    ism1 = akli @ ism2
    ism2 = ism1 @ lk
    ogl -= 2.0 * vech(ism2, ne)

    dot1_st1 = st1.reshape(ne, 3).T
    ism1 = dot2_ee @ dot1_st1
    ism2 = akli @ ism1
    ism3 = ism2 @ lk
    ogl += 2.0 * vech(ism3, ne)

    ism3_T = ism2.T
    ism1 = ism3_T @ lk
    ogl += 2.0 * vech(ism1, ne)

    ogl *= ov

    # ------------------------------------------------------------------
    # Kinetic-energy integral and its gradients
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

    # Gradient of T w.r.t. sk  (kgs)
    al_mass = sh.mass @ al
    aklp = akb @ (massb @ alb)
    aklpt = aklp.T

    im1_tmp = aklib @ akb
    im2_tmp = aklp @ im1_tmp
    im2_tmp = im2_tmp.T
    iv2_kgs = im2_tmp @ dif1

    im1_tmp = aklib @ akb
    im2_tmp = aklpt @ im1_tmp
    im3_tmp = im2_tmp - aklpt
    im3_tmp = im3_tmp.T
    iv1_kgs = im3_tmp @ dif2

    kgs = (2.0 * (2.0 * ke_i3 + ke_i2) * ogs
           + 4.0 * ov * iv2_kgs
           + 4.0 * ov * iv1_kgs)

    # Gradient of T w.r.t. vech(lk)  (kgl)
    kgl = (2.0 * ke_i3 + ke_i2) * ogl

    ism1 = akli @ lk
    ism2 = al_mass @ ism1
    ism3 = ak @ ism2
    ism1 = akli @ ism3
    kgl -= ov * 6.0 * vech(ism1, ne)

    ism1 = akli @ lk
    ism2 = al_mass @ ism1
    kgl += ov * 6.0 * vech(ism2, ne)

    # Shape matrices for the s, sk, sl quantities
    dot1_s = st.reshape(ne, 3).T
    dot2_s = dot1_s.T
    dot1_sk = st1.reshape(ne, 3).T
    dot2_sk = dot1_sk.T

    sds = dot2_s @ dot1_s
    sdsk = dot2_s @ dot1_sk
    skdsk = dot2_sk @ dot1_sk

    dot1_sl = st2.reshape(ne, 3).T
    dot2_sl = dot1_sl.T
    sldsk = dot2_sl @ dot1_sk
    slds = dot2_sl @ dot1_s

    aklps = np.kron(ak @ al_mass, I3)
    aklpst = aklps.T

    # kcm accumulation (ne × ne)
    def _ne_mat(m_kron: np.ndarray) -> np.ndarray:
        """Extract ne×ne matrix from block-diagonal skp×skp matrix."""
        # Undo kron(A, I3): A[i,j] = m_kron[3*i, 3*j]
        out = np.zeros((ne, ne))
        for ii_ in range(ne):
            for jj_ in range(ne):
                out[ii_, jj_] = m_kron[3 * ii_, 3 * jj_]
        return out

    # The kcm terms are built in the (ne×ne) space directly to match
    # the Fortran loops which operate on ne×ne matrices before kron.
    al_ne = al_mass        # (ne, ne)
    aklp_ne = ak @ al_ne  # (ne, ne)  =  aklps before kron
    aklpt_ne = aklp_ne.T

    kcm = np.zeros((ne, ne))

    # term 1: (al - akli·(aklp+aklpt)) · sds
    t = al_ne - akli @ (aklp_ne + aklpt_ne)
    kcm += t @ sds

    # term 2: (akli·(aklp+aklpt) - al) · sdsk
    t = akli @ (aklp_ne + aklpt_ne) - al_ne
    kcm += t @ sdsk

    # term 3: akli·aklpt · sdsk^T
    sdsk_T = sdsk.T
    kcm += (akli @ aklpt_ne) @ sdsk_T

    # term 4: (akli·aklp - al) · slds
    kcm += (akli @ aklp_ne - al_ne) @ slds

    # term 5: (al - akli·aklp) · sldsk
    kcm += (al_ne - akli @ aklp_ne) @ sldsk

    # term 6: -akli·aklpt · skdsk
    kcm -= (akli @ aklpt_ne) @ skdsk

    ism1 = kcm @ lk
    kgl += 2.0 * ov * vech(ism1, ne)
    kcm_T = kcm.T
    ism1 = kcm_T @ lk
    kgl += 2.0 * ov * vech(ism1, ne)
    kgl *= 2.0

    # ------------------------------------------------------------------
    # Electron-repulsion gradients
    # ------------------------------------------------------------------
    ert = 0.0
    ergst = np.zeros(skp)
    erglt = np.zeros(svh)

    for pp in range(ne + 1):
        for ll_ in range(pp + 1, ne + 1):
            jij_orig = np.zeros((ne, ne))
            if pp == 0:
                jij_orig[ll_ - 1, ll_ - 1] = 1.0
            else:
                jij_orig[pp - 1, pp - 1] = 1.0
                jij_orig[ll_ - 1, ll_ - 1] = 1.0
                jij_orig[pp - 1, ll_ - 1] = -1.0
                jij_orig[ll_ - 1, pp - 1] = -1.0

            ism1 = jij_orig @ akli
            jij_new = akli @ ism1           # akli · jij_orig · akli
            jijb = np.kron(jij_new, I3)

            iv3_er = jijb @ ee
            i3_er = float(ee @ iv3_er)
            i2_er = float(np.trace(ism1))  # tr(jij_orig · akli)
            if i2_er == 0.0:
                continue
            i1_er = i3_er / i2_er

            er = (sh.ncv[pp] * sh.ncv[ll_]
                  * 2.0 * ov * f0(i1_er) / math.sqrt(i2_er * sh.PI))
            ert += er

            # Gradient of er w.r.t. sk
            jijb_tmp = akb @ jijb   # used for iv3 recompute
            iv3_tmp = jijb_tmp @ ee
            ergs = ((sh.ncv[pp] * sh.ncv[ll_] * 2.0 * f0(i1_er) / math.sqrt(i2_er * sh.PI)) * ogs
                    + sh.ncv[pp] * sh.ncv[ll_] * 4.0 * ov * df0(i1_er) * iv3_tmp / math.sqrt(sh.PI * i2_er ** 3))
            ergst += ergs

            # Gradient of er w.r.t. vech(lk)
            ism1_jij = jij_new @ lk    # akli · jij_orig · akli · lk
            ergl = ((2.0 * f0(i1_er) / math.sqrt(i2_er * sh.PI)) * ogl
                    + (2.0 * ov * f0(i1_er) / math.sqrt(sh.PI * i2_er ** 3)) * vech(ism1_jij, ne))

            ergl += 4.0 * ov * df0(i1_er) * i3_er * vech(ism1_jij, ne) / math.sqrt(sh.PI * i2_er ** 5)

            # m2 term
            dot1_ee2 = ee.reshape(ne, 3).T
            dot2_ee2 = dot1_ee2.T
            ism1_m2 = dot2_ee2 @ dot1_ee2
            ism2_m2 = akli @ ism1_m2
            ism1_m2 = ism2_m2 @ jij_new

            ism2_m2 = ism1_m2 @ lk
            ergl -= 4.0 * ov * df0(i1_er) * vech(ism2_m2, ne) / math.sqrt(sh.PI * i2_er ** 3)

            ism2_m2T = ism1_m2.T
            ism1_m2 = ism2_m2T @ lk
            ergl -= 4.0 * ov * df0(i1_er) * vech(ism1_m2, ne) / math.sqrt(sh.PI * i2_er ** 3)

            # m1 term
            dot1_sk2 = st1.reshape(ne, 3).T
            dot2_sk2 = dot1_sk2.T
            ism1_m1 = dot2_sk2 @ dot1_ee2
            ism2_m1 = ism1_m1 @ jij_new

            ism1_m1 = ism2_m1 @ lk
            ergl += 4.0 * ov * df0(i1_er) * vech(ism1_m1, ne) / math.sqrt(sh.PI * i2_er ** 3)

            ism1_m1T = ism2_m1.T
            ism2_m1 = ism1_m1T @ lk
            ergl += 4.0 * ov * df0(i1_er) * vech(ism2_m1, ne) / math.sqrt(sh.PI * i2_er ** 3)

            erglt += sh.ncv[pp] * sh.ncv[ll_] * ergl

    # ------------------------------------------------------------------
    # Electric-field (Stark) gradient
    # ------------------------------------------------------------------
    zht_sum = 0.0
    for i in range(1, ne + 1):
        zht_sum += sh.ncv[i] * st[3 * i - 1]

    zgs = np.zeros(skp)
    im1_z = akb @ aklib
    for i in range(1, ne + 1):
        zgs += sh.ncv[i] * (ogs * st[3 * i - 1] + ov * im1_z[:, 3 * i - 1]) * sh.efs

    zgl = np.zeros(svh)
    voo = np.zeros(skp)
    for kk in range(1, ne + 1):
        voo[:] = 0.0
        voo[3 * kk - 1] = 1.0   # z-component unit vector for particle kk
        dot1_z = (st1 - st).reshape(ne, 3).T
        dot1e_z = voo.reshape(ne, 3).T

        dot2_z = dot1_z.T
        ism1_z = dot2_z @ dot1e_z
        ism2_z = ism1_z @ akli
        ism2_z = ism2_z @ lk
        v1_z = vech(ism2_z, ne)

        ism1_zT = ism1_z.T
        ism2_zT = akli @ ism1_zT
        ism2_zT = ism2_zT @ lk
        v2_z = vech(ism2_zT, ne)

        zgl += sh.ncv[kk] * sh.efs * (st[3 * kk - 1] * ogl + ov * (v1_z + v2_z))

    # ------------------------------------------------------------------
    # Total gradient vectors
    # igvl = kgl + erglt - zgl
    # igvs = kgs + ergst - zgs
    # ------------------------------------------------------------------
    igvl = kgl + erglt - zgl
    igvs = kgs + ergst - zgs

    tgvl = np.empty(svh + skp)
    tgvls = np.empty(svh + skp)
    tgvl[:svh] = igvl
    tgvl[svh:] = igvs
    tgvls[:svh] = ogl
    tgvls[svh:] = ogs

    # Diagonal elements get doubled (bra == ket)
    if ii == jj:
        tgvl *= 2.0
        tgvls *= 2.0

    return tgvl, tgvls
