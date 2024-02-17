#ifndef CHOLESKY_H
#define CHOLESKY_H
#include <vector>
#include <math.h>

#include <spamtrix_setup.hpp>
#include <spamtrix_ircmatrix.hpp>
#include <spamtrix_fleximatrix.hpp>


namespace SpaMtrix {
class Cholesky {
  FlexiMatrix L;    // Lower diagonal matrix

  Cholesky():L(){}
  void forwardSubstitution(Vector&x, const Vector& b) const;
  void backwardSubstitution(Vector&x, const Vector& b) const;

public:
    /** Creates a cholesky decomposition lower diagonal of matrix A*/
    explicit Cholesky(const IRCMatrix& A);
    void print()const;

    /**
        Solves Ax=b using forward/backward substitution. <br>
        Ax = b <br>
        L'(Ly) = b <br>
        L'x = y <br>
    */
    void solve(Vector& x, const Vector& b) const; // SOLVES Ax=b USING FORWARD/BACKWARD SUBSTITUTION
    virtual ~Cholesky();
};
} // end namespace SpaMtrix

#endif // CHOLESKY_H
