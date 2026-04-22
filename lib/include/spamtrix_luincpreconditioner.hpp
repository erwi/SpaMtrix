// Copyright (c) 2026 Eero
// SPDX-License-Identifier: MIT
#ifndef LUINCPRECONDITIONER_H
#define LUINCPRECONDITIONER_H

#include <spamtrix_preconditioner.hpp>
#include <spamtrix_fleximatrix.hpp>
#include <spamtrix_ircmatrix.hpp>
#include <spamtrix_vector.hpp>

namespace SpaMtrix {
    /**
     * @brief Incomplete LU preconditioner.
     */
    class LUIncPreconditioner: public Preconditioner {
        FlexiMatrix M;  ///< Stores both L and U factors.
        LUIncPreconditioner(){}

        /**
            @brief Forward substitution on the lower triangular factor.
            
            The lower factor is unit diagonal, so the diagonal is not stored.
            
            @param x Solution vector to populate.
            @param b Right-hand-side vector.
        */
        void forwardSubstitution(Vector &x, const Vector &b) const;

        /**
            @brief Back-substitution on the upper triangular factor.
            
            @param x Solution vector to populate.
            @param b Right-hand-side vector.
        */
        void backwardSubstitution(Vector &x, const Vector &b) const;
    public:
        /**
         * @brief Create an incomplete LU preconditioner with zero drop tolerance.
         *
         * @param A Source matrix.
         */
        LUIncPreconditioner(const IRCMatrix &A);
        /** @brief Destroy the preconditioner. */
        virtual ~LUIncPreconditioner();
        /** @brief Print the factorization. */
        void print() const;
        /**
         * @brief Solve the ILU preconditioner system.
         *
         * @param x Solution vector to populate.
         * @param b Right-hand-side vector.
         */
        void solveMxb(Vector &x, const Vector &b) const;
    };
}
#endif
