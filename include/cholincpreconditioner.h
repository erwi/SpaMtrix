#ifndef CHOLINCPRECONDITIONER_H
#define CHOLINCPRECONDITIONER_H
#include <setup.h>
#include <ircmatrix.h>
#include <vector.h>
#include <fleximatrix.h>
#include <preconditioner.h>

#include <iostream>
#include <math.h>

namespace SpaMtrix
{

class CholIncPreconditioner:public Preconditioner
{
  FlexiMatrix L;

  CholIncPreconditioner():L(0,0){}
  void forwardSubstitution(Vector&x, const Vector& b) const;
  void backwardSubstitution(Vector&x, const Vector& b) const;
public:
  CholIncPreconditioner(const IRCMatrix &A);
  void print() const;
  void solveMxb(Vector &x, const Vector &b) const;
};
} // end namespace SpaMtrix
#endif
