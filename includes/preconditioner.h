#ifndef PRECONDITIONER_H
#define PRECONDITIONER_H
#include <vector.h>
// PURE VIRTUAL PRECONDITIONER BASE CLASS
class Preconditioner
{
public:
  virtual void solveMxb(Vector &x, const Vector &b) const  = 0; // SOLVES Mx = b
};

#endif
