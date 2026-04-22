// Copyright (c) 2026 Eero
// SPDX-License-Identifier: MIT
#ifndef DIAGPRECONDITIONER_H
#define DIAGPRECONDITIONER_H
#include <omp.h>

#include <spamtrix_setup.hpp>
#include <spamtrix_vector.hpp>
#include <spamtrix_ircmatrix.hpp>
#include <spamtrix_preconditioner.hpp>


namespace SpaMtrix
{
/**
 * @brief Diagonal preconditioner built from the diagonal of an IRC matrix.
 */
class DiagPreconditioner: public Preconditioner {

    Vector diagonal;

    DiagPreconditioner();
public:
    /**
     * @brief Construct the diagonal preconditioner from a sparse matrix.
     *
     * @param A Source matrix.
     */
    DiagPreconditioner(const IRCMatrix& A);
    /** @brief Destroy the preconditioner. */
    virtual ~DiagPreconditioner();
    /**
     * @brief Solve the diagonal preconditioner system.
     *
     * @param x Solution vector to populate.
     * @param b Right-hand-side vector.
     */
    void solveMxb(Vector &x, const Vector &b) const;
};
} // end namespace SpaMtrix

#endif

