#ifndef JACOBIPRECONDITIONER_H
#define JACOBIPRECONDITIONER_H
#include <setup.h>
#include <preconditioner.h>
#include <ircmatrix.h>
#include <vector.h>
#include <relaxers.h>
class JacobiPreconditioner:public Preconditioner
{
  IRCMatrix const * const A; 	// POINTER TO MATRIX THAT IS PRECONDITIONED
  const idx numIterations;	// NUMBER OF JACOBI ITERATIONS PERFORMED

public:
  JacobiPreconditioner(const IRCMatrix& A, const idx numIterations = 1):
  A( &A),
  numIterations(numIterations)
  {  }
  
  void solveMxb(Vector &x, const Vector &b) const
  {
    jacobi(*A , x , b, numIterations );
  }
  
};

#endif

