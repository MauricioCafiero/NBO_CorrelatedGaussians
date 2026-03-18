# NECG — N-electron Correlated Gaussians

A quantum-chemistry program for computing energies and properties of small
atomic and molecular systems using **explicitly-correlated Gaussian (ECG)**
basis functions in the **non-Born-Oppenheimer** framework.

Originally written in Fortran with MPI parallelism (see `source/`), and
subsequently refactored as a single-process Python package (see `python/`).

---

## Physics background

The wave function is expressed as a linear combination of correlated Gaussians
with shifted centres:

```math
\Psi = \sum_k C_k \, \varphi_k
```

```math
\varphi_k = \exp\!\bigl(-\mathbf{r}^T A_k \mathbf{r}\bigr)
            \exp\!\bigl(-\mathbf{s}_k^T \mathbf{r}\bigr)
```

where $A_k = L_k L_k^T$ is a positive-semidefinite matrix constructed from
a lower-triangular Cholesky factor $L_k$, and $\mathbf{s}_k$ is a shift
vector.  Both the nonlinear parameters ($L_k$, $\mathbf{s}_k$) and the
linear expansion coefficients ($C_k$) are optimised by minimising the
Rayleigh quotient

```math
E = \frac{\langle C | H | C \rangle}{\langle C | S | C \rangle}
```

where $H$ is the Hamiltonian matrix and $S$ is the overlap matrix.

The approach captures inter-particle correlation explicitly and does not
invoke the Born-Oppenheimer approximation, making it well-suited to
high-accuracy calculations on small systems such as He, $H_2$, and $Li_2$.

---

## Repository layout

```
NBO_CorrelatedGaussians/
├── source/          # Original Fortran 90 / Fortran 77 code (MPI)
├── python/          # Python refactoring (single-process, NumPy/SciPy)
├── oldquit.dat      # Default H2 input file (16-term basis)
├── h2_input         # Alternative H2 example input
└── LICENSE          # MIT License
```

---

## Original Fortran code (`source/`)

### Overview

The Fortran implementation is the original production code, designed for
parallel execution on HPC clusters using MPI.  It is compiled with
`mpif77` and links against LAPACK and BLAS.

**Attribution:** Adamowicz Research Group, University of Arizona Chemistry.
Programming by Mauricio Cafiero.

### Building

```bash
cd source
make
```

This produces the `necg` executable.  The `makefile` requires:

* An MPI Fortran compiler (`mpif77`)
* LAPACK and BLAS libraries (`-llapack -lblas`)

### Running

```bash
mpirun -np <nproc> ./necg
```

The program reads `oldquit.dat` (or the file path specified in the source)
and writes results to `necg.out`.

### Source files

| File | Purpose |
|------|---------|
| `necg.f90` | Main driver; initialises matrices, calls optimiser |
| `matel.f90` | Matrix elements (overlap, kinetic energy, electron repulsion) |
| `gradel.f90` | Analytical gradients of all matrix elements |
| `make_grad.f90` | Assembles the full gradient vector |
| `outprod.f90` | Rayleigh-quotient objective function; calls gradient assembly |
| `rrms.f90` | Dipole-moment matrix elements |
| `read_input.f90` | Reads the `oldquit.dat` input file |
| `matrices_to_be_shared.f90` | Global shared-state module (MPI) |
| `tn.f` | Truncated-Newton (LMQN) optimiser |
| `btn.f` | Alternative BFGS optimiser |
| `inv.f90` | Matrix inversion via LAPACK `DPOTRF`/`DPOTRI` |
| `vech.f90` | Vector/half-vectorisation utilities |
| `biout.f90` | Formatted output helper |
| `matrix_print.f90` | Matrix printing utility |
| `error_function_approximation.f90` | Padé approximation for the Boys function $F_0(z)$ |
| `derror.f90` | Derivative of the Boys function |
| `get_current_date_and_time.f90` | Timestamp utility |

### Key computational features

* **MPI parallelism** — matrix-element and gradient loops are distributed
  across multiple MPI ranks.
* **LAPACK/BLAS** — `DPOTRF`/`DPOTRI` for Cholesky inversion; `DSYGV` for
  the generalised eigenvalue problem.
* **Boys function** — electron-repulsion integrals use the Padé approximation
  for $F_0(z) = \int_0^1 e^{-z t^2} dt$.
* **Cholesky parameterisation** — guarantees positive-semidefiniteness of
  $A_k$ without constraints.

---

## Python refactoring (`python/`)

### Overview

A direct, module-for-module port of the Fortran code to Python.  MPI
parallelism has been removed; all loops run sequentially on a single
process.  The input/output format is identical to the Fortran version.

See [`python/README.md`](python/README.md) for full details.

### Requirements

```bash
pip install -r python/requirements.txt
```

Dependencies: **NumPy** ≥ 1.24, **SciPy** ≥ 1.10.

### Running

```bash
cd python
python necg.py [input_file]
```

`input_file` defaults to `../oldquit.dat`.  Output is written to stdout
and to `necg.out`.

### Module map

| Python file | Fortran analogue | Purpose |
|-------------|-----------------|---------|
| `shared.py` | `matrices_to_be_shared.f90` | Global state (dimensions, matrices, vectors) |
| `utils.py` | `biout.f90`, `matrix_print.f90`, `error_function_approximation.f90`, `derror.f90`, `get_current_date_and_time.f90`, `vech.f90` | Utility functions |
| `read_input.py` | `read_input.f90` | Read `oldquit.dat` |
| `matel.py` | `matel.f90` | Matrix elements (overlap, kinetic, repulsion) |
| `gradel.py` | `gradel.f90` | Analytical gradients of matrix elements |
| `make_grad.py` | `make_grad.f90` | Assemble full gradient vector |
| `outprod.py` | `outprod.f90` | Objective function + gradient for optimiser |
| `rrms.py` | `rrms.f90` | Dipole-moment matrix elements |
| `necg.py` | `necg.f90` | Main driver |

### Key differences from Fortran

* **MPI removed** — all parallelism dropped; runs on a single process.
* **NumPy for linear algebra** — `DPOTRF`/`DPOTRI` → `numpy.linalg.inv`;
  `DSYGV` → `scipy.linalg.eigh`; Kronecker products →
  `numpy.kron(A, numpy.eye(3))`.
* **SciPy optimiser** — the Truncated-Newton / LMQN routine (`tn.f`) is
  replaced by `scipy.optimize.minimize(method="L-BFGS-B")`.
* **Boys function** — Padé approximation from the Fortran source preserved
  exactly.
* **Input/output format** unchanged — `oldquit.dat` is read with the same
  layout; `necg.out` has the same structure.

---

## Input / output format

### Input (`oldquit.dat`)

A fixed-layout text file.  The first few records control the basis and
system:

| Record | Description |
|--------|-------------|
| `M` | Number of Gaussian basis functions |
| `ne` | Number of electrons |
| `nst` | Number of symmetry terms |
| `efs` | Electric-field strength (0 = no field) |
| `mass(ne×ne)` | Electron mass tensor |
| `ncv(ne+1)` | Nuclear charges |
| `iv(M×(svh+skp))` | Nonlinear parameters (Cholesky factors + shift vectors) |
| `C(M)` | Linear expansion coefficients |
| `scvv(nst)` | Symmetry-term weights |
| `sml`, `sms` | Symmetry matrices for Gaussians and shifts |

Example input files for $H_2$ are provided as `oldquit.dat` (16-term basis)
and `h2_input`.

### Output (`necg.out`)

| Section | Contents |
|---------|---------|
| Header | Program title, date/time, parameter summary |
| Initial matrices | $H$, $S$, kinetic and potential energy matrices |
| Initial energy | Rayleigh quotient before optimisation |
| Optimisation log | L-BFGS-B (or LMQN) convergence iterations |
| Final energy | Optimised ground-state energy (hartree) |
| Final properties | Average kinetic/potential energies, virial coefficient, dipole moment (a.u. and Debye) |

---

## License

MIT — see [LICENSE](LICENSE).
