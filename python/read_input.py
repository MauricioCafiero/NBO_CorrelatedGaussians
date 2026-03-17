"""
Read the input file (oldquit.dat) and populate the shared global state.

Analogous to read_input.f90.  MPI broadcast calls have been removed; all
data is read on a single process.
"""

from __future__ import annotations

import numpy as np

import shared as sh
from utils import biout


def read_input(filename: str = "oldquit.dat") -> None:
    """
    Read *filename* and populate every field in *shared*.

    File layout (one value or row of values per line)::

        M
        ne
        nst
        efs
        mass(1,1) mass(1,2) ... mass(1,ne)
        mass(2,1) ...
        ...
        mass(ne,1) ...
        ncv(1)
        ncv(2)
        ...
        ncv(ne+1)
        iv(1)
        iv(2)
        ...
        iv(siv)
        C(1)
        ...
        C(M)
        scvv(1)
        ...
        scvv(nst)
        sml(1,1,1) sml(1,2,1) ... sml(1,ne,1)   ! row 1 of sml(:,:,1)
        sml(2,1,1) ...
        ...                                         ! all ne rows for k=1
        sml(1,1,2) ...                              ! k=2
        ...                                         ! repeat for all nst
        sms(1,1,1) sms(1,2,1) ... sms(1,skp,1)
        sms(2,1,1) ...
        ...                                         ! all skp rows for k=1
        ...                                         ! repeat for all nst
    """
    try:
        with open(filename, "r", encoding="utf-8") as fh:
            tokens = fh.read().split()
    except OSError:
        biout(f"Could not open input file: {filename}")
        raise

    pos = 0

    def _next_int() -> int:
        nonlocal pos
        val = int(tokens[pos])
        pos += 1
        return val

    def _next_float() -> float:
        nonlocal pos
        val = float(tokens[pos].lower().replace('d', 'e'))
        pos += 1
        return val

    # ------------------------------------------------------------------
    # Scalar parameters
    # ------------------------------------------------------------------
    sh.M = _next_int()
    sh.ne = _next_int()
    sh.nst = _next_int()
    sh.efs = _next_float()

    sh.skp = 3 * sh.ne
    sh.svh = sh.ne * (sh.ne + 1) // 2
    sh.siv = sh.M * (sh.svh + sh.skp)

    # ------------------------------------------------------------------
    # Allocate all arrays
    # ------------------------------------------------------------------
    sh.mass = np.zeros((sh.ne, sh.ne))
    sh.ki = np.zeros((sh.M, sh.M))
    sh.po = np.zeros((sh.M, sh.M))
    sh.hh = np.zeros((sh.M, sh.M))
    sh.ss = np.zeros((sh.M, sh.M))
    sh.C = np.zeros(sh.M)
    sh.S = np.zeros((sh.M, sh.M))
    sh.H = np.zeros((sh.M, sh.M))
    sh.Kinetic = np.zeros((sh.M, sh.M))
    sh.Potential = np.zeros((sh.M, sh.M))
    sh.grad = np.zeros(sh.siv + sh.M)
    sh.DS = np.zeros((sh.M * (sh.svh + sh.skp), sh.M))
    sh.DH = np.zeros((sh.M * (sh.svh + sh.skp), sh.M))
    sh.DSS = np.zeros((sh.M * (sh.svh + sh.skp), sh.M))
    sh.DHH = np.zeros((sh.M * (sh.svh + sh.skp), sh.M))
    sh.ncv = np.zeros(sh.ne + 1)
    sh.scvv = np.zeros(sh.nst)
    sh.iv = np.zeros(sh.siv)
    sh.sml = np.zeros((sh.ne, sh.ne, sh.nst))
    sh.sms = np.zeros((sh.skp, sh.skp, sh.nst))

    # ------------------------------------------------------------------
    # Mass matrix  (ne × ne, stored row-by-row in input file)
    # ------------------------------------------------------------------
    for i in range(sh.ne):
        for j in range(sh.ne):
            sh.mass[i, j] = _next_float()

    # ------------------------------------------------------------------
    # Charge values ncv(1..ne+1)
    # ------------------------------------------------------------------
    for i in range(sh.ne + 1):
        sh.ncv[i] = _next_float()

    # ------------------------------------------------------------------
    # Nonlinear parameter vector iv(1..siv)
    # ------------------------------------------------------------------
    for i in range(sh.siv):
        sh.iv[i] = _next_float()

    # ------------------------------------------------------------------
    # Linear coefficients C(1..M)
    # ------------------------------------------------------------------
    for i in range(sh.M):
        sh.C[i] = _next_float()

    # ------------------------------------------------------------------
    # Symmetry weights scvv(1..nst)  — normalised by nst
    # ------------------------------------------------------------------
    for i in range(sh.nst):
        sh.scvv[i] = _next_float()
    sh.scvv /= sh.nst

    # ------------------------------------------------------------------
    # Symmetry matrix sml(ne, ne, nst)  — stored row-by-row per k
    # ------------------------------------------------------------------
    for k in range(sh.nst):
        for i in range(sh.ne):
            for j in range(sh.ne):
                sh.sml[i, j, k] = _next_float()

    # ------------------------------------------------------------------
    # Symmetry matrix for shifts sms(skp, skp, nst)
    # ------------------------------------------------------------------
    for k in range(sh.nst):
        for i in range(sh.skp):
            for j in range(sh.skp):
                sh.sms[i, j, k] = _next_float()
