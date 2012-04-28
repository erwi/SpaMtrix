#ifndef SORPRECONDITIONER_H
#define SORPRECONDITIONER_H
#include <setup.h>
#include <ircmatrix.h>
#include <vector.h>
#include <preconditioner.h>

class SORPreconditioner: public Preconditioner
{
  IRCMatrix const * const A;
  const idx numIterations;
  
public:
  SORPreconditioner(const IRCMatrix& A, const idx numIterations = 1):
  A(&A),
  numIterations(numIterations){}
  
  void solveMxb(Vector &x, const Vector &b) const
  {
    SOR(*A , x, b, numIterations );
  }
  
};
#endif
