#ifndef CHOLESKY_H
#define CHOLESKY_H
#include <vector>
#include <math.h>

#include <spamtrix_setup.hpp>
#include <spamtrix_ircmatrix.hpp>
#include <spamtrix_fleximatrix.hpp>


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
    virtual ~Cholesky();
};
} // end namespace SpaMtrix

#endif // CHOLESKY_H
