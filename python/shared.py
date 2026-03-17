"""
Shared global state for N-electron Correlated Gaussians (NECG) calculations.

Analogous to matrices_to_be_shared.f90.  All mutable state is stored as
module-level variables so that it can be imported and mutated by every
sub-module, mirroring the Fortran USE association mechanism.
"""

import numpy as np

# ---------------------------------------------------------------------------
# Dimensional parameters (set by read_input)
# ---------------------------------------------------------------------------
M: int = 0      # number of basis functions (size of Hamiltonian matrix)
ne: int = 0     # number of electrons
NN: int = 0     # number of nuclei
nst: int = 0    # number of symmetry terms
svh: int = 0    # ne*(ne+1)/2  (size of vech of Ak)
skp: int = 0    # 3*ne          (Kronecker-product dimension)
siv: int = 0    # M*(svh+skp)   (total size of input vector iv)

# ---------------------------------------------------------------------------
# Physical / field parameters
# ---------------------------------------------------------------------------
efs: float = 0.0   # electric-field strength

# ---------------------------------------------------------------------------
# Matrices allocated after M, ne are known
# ---------------------------------------------------------------------------
mass: np.ndarray = None       # (ne, ne)  inverse-mass tensor
S: np.ndarray = None          # (M, M)    overlap matrix
H: np.ndarray = None          # (M, M)    Hamiltonian matrix
Kinetic: np.ndarray = None    # (M, M)    kinetic-energy matrix
Potential: np.ndarray = None  # (M, M)    potential-energy matrix
ki: np.ndarray = None         # (M, M)    work copy (kinetic)
po: np.ndarray = None         # (M, M)    work copy (potential)
hh: np.ndarray = None         # (M, M)    work copy (Hamiltonian)
ss: np.ndarray = None         # (M, M)    work copy (overlap)
DS: np.ndarray = None         # (M*(svh+skp), M)  gradient of S
DH: np.ndarray = None         # (M*(svh+skp), M)  gradient of H
DSS: np.ndarray = None        # (M*(svh+skp), M)  reduced gradient of S
DHH: np.ndarray = None        # (M*(svh+skp), M)  reduced gradient of H
GRADM: np.ndarray = None      # (M, M*(svh+skp))  gradient matrix
grad: np.ndarray = None       # (M*(svh+skp+1),)  full gradient vector
GRADA: np.ndarray = None      # (M*(svh+skp),)    nonlinear-param gradient
GRADC: np.ndarray = None      # (M,)              linear-coeff gradient

# ---------------------------------------------------------------------------
# Vectors allocated after M, ne are known
# ---------------------------------------------------------------------------
iv: np.ndarray = None         # (siv,)   nonlinear parameters (Cholesky + shifts)
ncv: np.ndarray = None        # (ne+1,)  charge values
scvv: np.ndarray = None       # (nst,)   symmetry weights
C: np.ndarray = None          # (M,)     linear coefficients

# ---------------------------------------------------------------------------
# 3-D symmetry arrays
# ---------------------------------------------------------------------------
sml: np.ndarray = None        # (ne, ne, nst)    symmetry matrix for Ak
sms: np.ndarray = None        # (skp, skp, nst)  symmetry matrix for shifts

# ---------------------------------------------------------------------------
# Scalar used by LAPACK-equivalent routines (Cholesky info flag)
# ---------------------------------------------------------------------------
jjj: int = 0

# ---------------------------------------------------------------------------
# Numerical constants
# ---------------------------------------------------------------------------
PI: float = 3.1415926535898

ZERO: float = 0.00
HALF: float = 0.50
QUARTER: float = 0.25
EIGHTH: float = 0.125
ONE: float = 1.00
TWO: float = 2.00
THREE: float = 3.00
FOUR: float = 4.00
FIVE: float = 5.00
SIX: float = 6.00
SEVEN: float = 7.00
EIGHT: float = 8.00
SIXTEEN: float = 16.00

# ---------------------------------------------------------------------------
# Pade-approximation constants for the Boys function F0
# ---------------------------------------------------------------------------
C1: float = 0.564691197e-04
C2: float = 0.758433197e-03
C3: float = 0.769838037e-02
C4: float = 0.629344460e-01
C5: float = 0.213271302e+00
C6: float = 0.720266520e-04
C7: float = 0.955528842e-03
C8: float = 0.101431553e-01
C9: float = 0.738522953e-01
C10: float = 0.338450368e+00
C11: float = 0.879937801e+00

# ---------------------------------------------------------------------------
# Fixed nuclear positions (hard-coded in the original)
# ---------------------------------------------------------------------------
EZ1: float = 0.0
EX1: float = 0.0
EY1: float = 0.952627944162921

EZ2: float = 0.0
EZ3: float = 0.0
EX2: float = 0.825000000000027
EX3: float = -0.825000000000027
EY2: float = -0.476313972081468
EY3: float = -0.476313972081468
