# Input File Generation for LiH

This directory contains tools for generating NECG input files for the LiH
(Lithium Hydride) molecule with any desired number of correlated Gaussian
basis functions.

---

## Background

The reference file `oldquit.dat` (in the repository root) encodes an optimised
12-function explicitly-correlated Gaussian basis for LiH in the non-Born-Oppenheimer framework.  For each basis function the file stores:

* A lower-triangular Cholesky factor **L** (packed as a vector, 15 elements
  for the 5-particle LiH system) such that **A** = **L L**ᵀ is the positive-
  definite exponent matrix of the Gaussian.
* A shift vector **s** (15 elements: x, y, z for each of the 5 pseudo-
  particles).

The generator reads these parameters from the reference file, computes their
per-column statistics, and samples novel but statistically similar parameters
for an arbitrary number of basis functions.

---

## Contents

| File | Purpose |
|------|---------|
| `generate_lih_input.py` | Main generator script |
| `README.md` | This file |

---

## Requirements

```bash
pip install numpy
```

NumPy is the only dependency (already listed in `python/requirements.txt`).

---

## Usage

### Command line

```bash
# Generate a 16-function LiH basis set (output: lih_M16.dat)
python generate_lih_input.py --M 16

# Specify the output filename explicitly
python generate_lih_input.py --M 8 --output lih_8.dat

# Fix the random seed for reproducibility
python generate_lih_input.py --M 20 --seed 42 --output lih_20.dat

# Use a custom reference file
python generate_lih_input.py --M 12 --reference /path/to/oldquit.dat

# Run the built-in self-tests (no --M required)
python generate_lih_input.py --test
```

### Python API

```python
from generate_lih_input import generate_lih_input

# Generate and write to a file; also returns the text
text = generate_lih_input(M=16, seed=42, output="lih_16.dat")

# Generate without writing (just get the string)
text = generate_lih_input(M=8, seed=0)
```

---

## Generated File Format

The output follows the exact same layout as `oldquit.dat` and can be used
directly as input to the NECG code:

```bash
# From the python/ directory:
python necg.py ../input_file_generation/lih_M16.dat
```

The fixed parameters (mass tensor, charges, symmetry matrices) are taken
unchanged from the reference file.  Only the basis-function–specific
nonlinear parameters (Cholesky elements + shift vectors) and linear
coefficients are newly generated.

---

## Generation Algorithm

1. **Parse** `oldquit.dat` to extract physical parameters and the 12×30
   matrix of nonlinear parameters (12 basis functions × 30 parameters each).

2. **Compute statistics**: per-column mean and standard deviation across the
   12 reference basis functions.

3. **Sample M new blocks**:
   * *Off-diagonal Cholesky elements* and *shift-vector components*: drawn
     from `Normal(mean, std)` for each column independently.
   * *Diagonal Cholesky elements* (must be strictly positive to maintain
     positive-definiteness): drawn from a `LogNormal` distribution whose
     mean and standard deviation match the reference column statistics.

4. **Set linear coefficients** to the equal-weight value `1/M`.

5. **Copy** all other parameters (mass tensor, charges, symmetry weights,
   and symmetry matrices) directly from the reference file.

The generated basis set is *novel* (sampled values differ from the reference)
but *statistically similar* (same distributional properties).  The NECG
optimiser will then refine all parameters to minimise the Rayleigh quotient.

---

## Self-tests

A built-in self-test suite verifies that:

* Generated files round-trip correctly through the parser.
* Header dimensions match the requested `M`.
* All diagonal Cholesky elements are strictly positive.
* Symmetry matrices and charge values are reproduced exactly.

```bash
python generate_lih_input.py --test
```

Expected output:
```
Running self-tests …
  PASS  M=  4: shape OK, diagonals positive, symmetry matrices exact
  PASS  M=  8: shape OK, diagonals positive, symmetry matrices exact
  PASS  M= 16: shape OK, diagonals positive, symmetry matrices exact
All self-tests passed.
```
