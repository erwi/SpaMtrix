# SpaMtrix
This is a small standalone library (i.e. doesn't depend on loads of, or even one, other libraries that must be first built and installed) for solving systems of linear simultaneous equations defined
using sparse matrices. The algorithms and datastructures are adapted from various books and other resources (e.g. Wikipedia) that may describe 
the method rather than best implementation, so are likely far from optimised.  

## Includes

Contains a sparse matrix class and a vector class for defining the system of linear equations Ax=b.

The following solver algorithms for solving the Ax=b system of equation are also implemented:

- Preconditioned Conjugate Gradient (PCG)
- Generalized Minimal Residual Method (GMRES)
- Cholesky decomposition method
- LU decomposition

The following preconditioners that can be used with the iterative solvers PCG and GMRES are implemented:

- Diagonal preconditioner (a.k.a. Jacobi preconditioner)
- Incomplete Cholesky factorisation preconditioner
- Incomplete LU factorisation preconditioner

## How to use it
The general idea is that first the matrix sparsity pattern is defined using the `MatrixMaker` class, which is a builder. This will then be converted
to the sparse matrix class SpaMtrix which can be used for actual computations. 

In absense of proper documentation, see the unit tests or examples in the `tests` or `examples` directories for details.

## Building 
Everything is built using CMake either as a library that can be linked to, or the whole projcet code can be included in your own CMake project. 
TODO: this can surely be improved.

This has been known to work at least with GCC and MinGW compilers, but there shouldn't be anything fundamentally preventing from using other complilers.



