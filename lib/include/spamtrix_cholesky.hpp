#ifndef CHOLESKY_H
#define CHOLESKY_H
#include <vector>
#include <math.h>

#include <spamtrix_setup.hpp>
#include <spamtrix_ircmatrix.hpp>
#include <spamtrix_fleximatrix.hpp>


namespace SpaMtrix {
/**
 * @brief Sparse Cholesky factorization.
 */
class Cholesky {
  FlexiMatrix L;

  Cholesky():L(){}
  /**
   * @brief Forward substitution with the lower factor.
   *
   * @param x Solution vector to populate.
   * @param b Right-hand-side vector.
   */
  void forwardSubstitution(Vector&x, const Vector& b) const;
  /**
   * @brief Backward substitution with the transposed lower factor.
   *
   * @param x Solution vector to populate.
   * @param b Right-hand-side vector.
   */
  void backwardSubstitution(Vector&x, const Vector& b) const;

public:
  /**
   * @brief Compute the Cholesky factorization of a sparse matrix.
   *
   * @param A Source matrix.
   */
    explicit Cholesky(const IRCMatrix& A);
  /** @brief Print the factorization. */
    void print()const;

    /**
   * @brief Solve $A x = b$ using forward and backward substitution.
   *
   * @param x Solution vector to populate.
   * @param b Right-hand-side vector.
    */
  void solve(Vector& x, const Vector& b) const;
  /** @brief Destroy the factorization. */
    virtual ~Cholesky();
};
} // end namespace SpaMtrix

#endif // CHOLESKY_H
