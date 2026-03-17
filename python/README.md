# NECG — N-electron Correlated Gaussians (Python)

Python refactoring of the Fortran source code located in `../source/`.

## Overview

The program computes energies and properties of small atomic / molecular
systems using explicitly-correlated Gaussian (ECG) basis functions in the
non-Born-Oppenheimer framework.  The wave function is a linear combination
of correlated Gaussians with shifted centres:

```
Ψ = Σ_k C_k · φ_k
φ_k = exp( -r^T A_k r ) · exp( -s_k^T r )
```

where `A_k = L_k L_k^T` is a positive-semidefinite matrix formed from a
lower-triangular Cholesky factor `L_k`.

## Requirements

```
pip install -r requirements.txt
```

Dependencies: **NumPy** ≥ 1.24, **SciPy** ≥ 1.10.

## Running

```
python necg.py [input_file]
```

The input file defaults to `oldquit.dat` (same format as the original
Fortran code).  Output is written to stdout and to `necg.out`.

## Module map

| Python file      | Fortran analogue              | Purpose                                      |
|------------------|-------------------------------|----------------------------------------------|
| `shared.py`      | `matrices_to_be_shared.f90`   | Global state (dimensions, matrices, vectors) |
| `utils.py`       | `biout.f90`, `matrix_print.f90`, `error_function_approximation.f90`, `derror.f90`, `get_current_date_and_time.f90`, `vech.f90` | Utility functions |
| `read_input.py`  | `read_input.f90`              | Read `oldquit.dat`                           |
| `matel.py`       | `matel.f90`                   | Matrix elements (overlap, kinetic, repulsion)|
| `gradel.py`      | `gradel.f90`                  | Analytical gradients of matrix elements      |
| `make_grad.py`   | `make_grad.f90`               | Assemble full gradient vector                |
| `outprod.py`     | `outprod.f90`                 | Objective function + gradient for optimiser  |
| `rrms.py`        | `rrms.f90`                    | Dipole-moment matrix elements                |
| `necg.py`        | `necg.f90`                    | Main driver                                  |

## Key changes from Fortran

* **MPI removed** — all parallelism has been dropped; loops run
  sequentially on a single process.
* **NumPy for linear algebra** — `DPOTRF`/`DPOTRI` → `numpy.linalg.inv`;
  `DSYGV` → `scipy.linalg.eigh`; Kronecker products →
  `numpy.kron(A, numpy.eye(3))`.
* **SciPy optimiser** — the Truncated-Newton / LMQN routine (`tn.f`) is
  replaced by `scipy.optimize.minimize(method="L-BFGS-B")`, a
  limited-memory quasi-Newton method with comparable convergence
  properties.
* **Boys function** — Pade approximation from
  `error_function_approximation.f90` / `derror.f90` preserved exactly.
* **Input / output format** unchanged — `oldquit.dat` is read with the
  same layout; `necg.out` has the same structure.
