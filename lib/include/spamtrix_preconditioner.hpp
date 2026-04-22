// Copyright (c) 2026 Eero
// SPDX-License-Identifier: MIT
#ifndef PRECONDITIONER_H
#define PRECONDITIONER_H
#include <spamtrix_vector.hpp>


namespace SpaMtrix {
/**
 * @brief Abstract base class for linear-system preconditioners.
 */
class Preconditioner {
    public:
        /**
         * @brief Solve the preconditioner system $M x = b$.
         *
         * @param x Solution vector to populate.
         * @param b Right-hand-side vector.
         */
        virtual void solveMxb(Vector &x, const Vector &b) const = 0;

        /**
         * @brief Solve the preconditioner system and return the result.
         *
         * @param b Right-hand-side vector.
         * @return The computed solution vector.
         */
        Vector solve(const Vector &b) const {
            Vector x(b.getLength());
            solveMxb(x, b);
            return x;
        }
        /** @brief Virtual destructor. */
        virtual ~Preconditioner() { }
    };

} // end namespace SpaMtrix

#endif
