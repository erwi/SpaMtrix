#ifndef CHOLESKY_H
#define CHOLESKY_H
#include <setup.h>
#include <ircmatrix.h>
#include <vector>
#include <math.h>
#include <fleximatrix.h>


namespace SpaMtrix
{
class Cholesky
{

  FlexiMatrix L;    // LOWER DIAGONAL


  Cholesky():L(0,0){}
  void forwardSubstitution(Vector&x, const Vector& b) const;
  void backwardSubstitution(Vector&x, const Vector& b) const;

public:
  Cholesky(const IRCMatrix& A);
  void print()const;

  void solve(Vector& x, const Vector& b) const; // SOLVES Ax=b USING FORWARD/BACKWARD SUBSTITUTION

};
} // end namespace SpaMtrix

#endif // CHOLESKY_H
