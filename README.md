# SpaMtrix

SpaMtrix is a small, standalone C++ library for solving systems of linear equations $Ax = b$ defined using sparse matrices. It has no external dependencies beyond a C++11-capable compiler and CMake. The algorithms and data structures are adapted from standard references and are intended to be straightforward rather than highly optimised.

## Contents

### Matrix types

| Class | Description |
|---|---|
| `IRCMatrix` | General sparse matrix using an Interleaved Row Compressed (IRC) storage scheme. This is the primary matrix type used for assembly and solving. |
| `FlexiMatrix` | Flexible sparse matrix used internally during matrix construction. |
| `TDMatrix` | Specialised tridiagonal matrix with a direct Thomas algorithm solver. |
| `DenseMatrix` | Dense matrix, used internally by the GMRES solver. |

### Vector

`Vector` — a dense vector class supporting standard arithmetic operators, norm computation, and normalisation.

### Solvers

**Iterative solvers** (via `IterativeSolvers`):
- Preconditioned Conjugate Gradient (PCG)
- Generalised Minimal Residual Method (GMRES)

**Direct solvers**:
- Cholesky decomposition (`Cholesky`) — for symmetric positive-definite matrices
- LU decomposition (`LU`) — for general square matrices (Crout algorithm, no pivoting)
- Tridiagonal direct solver (`TDMatrix::solveAxb`) — for tridiagonal systems

### Preconditioners

The following preconditioners can be used with the iterative solvers:

- `DiagPreconditioner` — diagonal (Jacobi) preconditioner
- `CholIncPreconditioner` — incomplete Cholesky factorisation preconditioner
- `LUIncPreconditioner` — incomplete LU factorisation preconditioner

### Eigenvalue computation

`powerMethod` — computes the dominant eigenvalue and corresponding eigenvector of a matrix using power iteration.

### I/O

- `Reader::readMatrixMarket` — reads a sparse matrix from a Matrix Market file
- `Writer::writeMatrixMarket` — writes a sparse matrix to a Matrix Market file
- `Writer::writeCSV` — writes a vector to a CSV file

## How to use it

The general workflow is:

1. Define the matrix sparsity pattern using `MatrixMaker`.
2. Populate entries with `addNonZero`.
3. Convert to `IRCMatrix` with `getIRCMatrix()`.
4. Construct a right-hand-side `Vector` and solve.

### Quick example: PCG solver

```cpp
#include <spamtrix_matrixmaker.hpp>
#include <spamtrix_ircmatrix.hpp>
#include <spamtrix_vector.hpp>
#include <spamtrix_iterativesolvers.hpp>
#include <spamtrix_diagpreconditioner.hpp>

using namespace SpaMtrix;

// Build a 2x2 matrix  A = [[3, 2], [2, 6]]
MatrixMaker mm(2, 2);
mm.addNonZero(0, 0, 3);
mm.addNonZero(0, 1, 2);
mm.addNonZero(1, 0, 2);
mm.addNonZero(1, 1, 6);
IRCMatrix A = mm.getIRCMatrix();

// Right-hand side b = [2, -8]
Vector b(2); b[0] = 2; b[1] = -8;
Vector x(2);

// Solve with PCG
DiagPreconditioner M(A);
IterativeSolvers solver(100, 1e-7);
bool converged = solver.pcg(A, x, b, M);
```

For more complete examples, see the `examples/` directory, and for tests see the `tests/` directory.

### Note on numeric types

`real` is a `typedef` for `double` and `idx` is a `typedef` for `unsigned int`. These are defined in `spamtrix_setup.hpp`.

## Building

### Prerequisites

- CMake 3.0 or later
- A C++11-capable compiler (GCC and MinGW are known to work; other compilers should work too)
- OpenMP (optional)

### Build steps

```bash
mkdir build && cd build
cmake ..
cmake --build .
```

### Optional: OpenMP support

Cholesky and LU decompositions have optional OpenMP parallelism. Enable it at configure time:

```bash
cmake -DUSES_OPENMP=ON ..
```

### Using SpaMtrix as a CMake sub-project

The library can be included directly in another CMake project:

```cmake
add_subdirectory(SpaMtrix)
target_link_libraries(MyTarget SpaMtrix)
```

The `lib/include` directory is automatically added to the include search path via the `SpaMtrix` target's public include directories.

### Installation

```bash
cmake --install .
```

This copies the library binary to `bin/` and the public headers to `include/`.

## Running tests

Tests use the [Catch2](https://github.com/catchorg/Catch2) (v1, single-header) framework, located in `extern/catch/`. Each Catch test case is registered as a separate CTest test, so `ctest` prints the individual test names as it runs them.

They are built automatically when SpaMtrix is the top-level CMake project and can be run with:

```bash
cd build
ctest --output-on-failure
# or directly
./tests/SpaMtrixTests
```

If your build directory is named differently, pass it explicitly instead of using `cd build`, for example:

```bash
ctest --test-dir cmake-build-debug --output-on-failure
```

For a more verbose run that includes extra CTest output, use `ctest -V --test-dir cmake-build-debug`.



