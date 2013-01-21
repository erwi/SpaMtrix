/* THIS IS A SWIG INTERFACE FILE FOR CREATING WRAPPERS FOR PYTHON */
%module SpaMtrix

%{
#include "../include/setup.h"
#include "../include/vector.h"
#include "../include/ircmatrix.h"
#include "../include/matrixmaker.h"
#include "../include/fleximatrix.h"
#include "../include/densematrix.h"
#include "../include/tdmatrix.h"
#include "../include/spamtrix_blas.h"
#include "../include/preconditioner.h"
#include "../include/diagpreconditioner.h"
#include "../include/cholincpreconditioner.h"
#include "../include/luincpreconditioner.h"
#include "../include/cholesky.h"
#include "../include/lu.h"
#include "../include/writer.h"
#include "../include/iterativesolvers.h"
%}
typedef unsigned int idx;
typedef double real;
// NEED TO RENAME OVERLOADED SPARSE MATRIX MULTIPLICATIONS
%rename(IRCMatMul) multiply(const IRCMatrix &, const Vector &, Vector &);
%rename(TDMatMul) multiply(const TDMatrix &, const Vector &, Vector &);

// IGNORE SOME OPERATORS THAT ARE NOT NEEDED IN PYTHON
%ignore multiply_dot(const IRCMatrix &, const Vector &, Vector &);
%ignore DenseMatrix::operator();
%ignore Vector::operator [];

%include "../include/setup.h"
%include "../include/vector.h"
%include "../include/ircmatrix.h"
%include "../include/matrixmaker.h"
%include "../include/fleximatrix.h"
%include "../include/densematrix.h"
%include "../include/tdmatrix.h"
%include "../include/spamtrix_blas.h"
%include "../include/preconditioner.h"
%include "../include/diagpreconditioner.h"
%include "../include/cholincpreconditioner.h"
%include "../include/luincpreconditioner.h"
%include "../include/cholesky.h"
%include "../include/lu.h"
%include "../include/writer.h"
%include "../include/iterativesolvers.h"
