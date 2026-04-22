// Copyright (c) 2026 Eero
// SPDX-License-Identifier: MIT
#ifndef SPAMTRIX_BLAS_H
#define SPAMTRIX_BLAS_H


#include <spamtrix_setup.hpp>


namespace SpaMtrix
{
class IRCMatrix;
class TDMatrix;
class Vector;

/**
 * @brief Compute the matrix-vector product of an IRC matrix and a vector.
 *
 * The result is written to @p b as $b = A x$.
 *
 * @param A Sparse IRC matrix.
 * @param x Input vector.
 * @param b Output vector receiving the product.
 */
void multiply(const IRCMatrix& A,
              const Vector& x,
              Vector& b);
/**
 * @brief Compute the matrix-vector product of a tridiagonal matrix and a vector.
 *
 * The result is written to @p b as $b = A x$.
 *
 * @param A Tridiagonal matrix.
 * @param x Input vector.
 * @param b Output vector receiving the product.
 */
void multiply(const TDMatrix& A,
              const Vector& x,
              Vector& b);
/**
 * @brief Compute $b = A x$ and return the dot product $x \cdot b$.
 *
 * This helper is used by iterative solvers that need both the matrix-vector
 * product and the resulting quadratic form.
 *
 * @param A Sparse IRC matrix.
 * @param x Input vector.
 * @param b Output vector receiving the product.
 * @return The dot product $x \cdot b$.
 */
real multiply_dot(const IRCMatrix& A,
                  const Vector& x,
                  Vector& b);

/**
 * @brief Return a vector containing the element-wise absolute values of @p vin.
 *
 * @param vin Input vector.
 * @return A copy of @p vin with each entry replaced by its absolute value.
 */
SpaMtrix::Vector abs(const SpaMtrix::Vector &vin);

/**
 * @brief Scale a vector in place.
 *
 * Applies $v \leftarrow a v$.
 *
 * @param a Scaling factor.
 * @param v Vector to modify.
 */
void scale(const real a, Vector& v);

/**
 * @brief Compute the dot product of two vectors.
 *
 * @param v1 First vector.
 * @param v2 Second vector.
 * @return The scalar product $v_1 \cdot v_2$.
 */
real dot(const Vector& v1, const Vector& v2);

/**
 * @brief Perform the BLAS operation $y \leftarrow y + a x$.
 *
 * @param a Scale factor for @p x.
 * @param x Input vector.
 * @param y Vector updated in place.
 */
void axpy(const real a, const Vector& x, Vector& y);

/**
 * @brief Perform the BLAS operation $y \leftarrow a y + x$.
 *
 * @param a Scale factor for @p y.
 * @param y Vector updated in place.
 * @param x Input vector.
 */
void aypx(const real a, Vector& y, const Vector& x);

/**
 * @brief Compute the squared residual norm for $A x = b$.
 *
 * Returns $\|b - A x\|_2^2$.
 *
 * @param A Sparse IRC matrix.
 * @param x Candidate solution vector.
 * @param b Right-hand-side vector.
 * @return The squared residual norm.
 */
real errorNorm2(const IRCMatrix& A, const Vector& x, const Vector& b);

/**
 * @brief Compute the squared residual norm for $A x = b$.
 *
 * Returns $\|b - A x\|_2^2$.
 *
 * @param A Tridiagonal matrix.
 * @param x Candidate solution vector.
 * @param b Right-hand-side vector.
 * @return The squared residual norm.
 */
real errorNorm2(const TDMatrix& A, const Vector& x, const Vector& b);

/**
 * @brief Return the Euclidean norm of a vector.
 *
 * @param x Input vector.
 * @return The 2-norm of @p x.
 */
real norm(const Vector& x);
} // end namespace SpaMtrix

#endif

