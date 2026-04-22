// Copyright (c) 2026 Eero
// SPDX-License-Identifier: MIT
#ifndef CHOLINCPRECONDITIONER_H
#define CHOLINCPRECONDITIONER_H
#include <iostream>
#include <math.h>

#include <spamtrix_setup.hpp>
#include <spamtrix_ircmatrix.hpp>
#include <spamtrix_vector.hpp>
#include <spamtrix_fleximatrix.hpp>
#include <spamtrix_preconditioner.hpp>

namespace SpaMtrix {
/**
 * @brief Incomplete Cholesky preconditioner.
 */
class CholIncPreconditioner: public Preconditioner {
    FlexiMatrix L;
    CholIncPreconditioner(): L() {}
    /**
     * @brief Forward substitution on the lower triangular factor.
     *
     * @param x Solution vector to populate.
     * @param b Right-hand-side vector.
     */
    void forwardSubstitution(Vector&x, const Vector& b) const;
    /**
     * @brief Back-substitution on the transposed lower triangular factor.
     *
     * @param x Solution vector to populate.
     * @param b Right-hand-side vector.
     */
    void backwardSubstitution(Vector&x, const Vector& b) const;
public:
    /**
     * @brief Construct the incomplete Cholesky preconditioner from a matrix.
     *
     * @param A Source matrix.
     */
    CholIncPreconditioner(const IRCMatrix &A);
    /** @brief Destroy the preconditioner. */
    virtual ~CholIncPreconditioner();
    /** @brief Print the factorization. */
    void print() const;
    /**
     * @brief Solve the incomplete Cholesky preconditioner system.
     *
     * @param x Solution vector to populate.
     * @param b Right-hand-side vector.
     */
    void solveMxb(Vector &x, const Vector &b) const;
};
} // end namespace SpaMtrix
#endif
