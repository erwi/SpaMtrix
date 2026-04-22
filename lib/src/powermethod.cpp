// Copyright (c) 2026 Eero
// SPDX-License-Identifier: MIT
#include <spamtrix_setup.hpp>
#include <spamtrix_blas.hpp>
#include <spamtrix_powermethod.hpp>
#include <spamtrix_vector.hpp>
#include <spamtrix_ircmatrix.hpp>
namespace SpaMtrix {
idx powerMethod(const SpaMtrix::IRCMatrix &A,
                real &eigenValue,
                SpaMtrix::Vector &eigenVector,
                real &toler,
                unsigned int maxIter
               ) {
    // Ensure the initial guess is not a zero vector.
    real n = eigenVector.getNorm();
    if (n == 0.0) {
        eigenVector(0) = 1.0;
    }
    eigenVector.normalise();
    real dL = toler + 1.0; // DELTA EIGENVALUE
    //real Lo(0);          // PREVIOUS EIGENVALUE
    SpaMtrix::Vector q(eigenVector);
    SpaMtrix::Vector z = A * eigenVector;
    idx iter(0);
    // Continue until the relative change falls below the tolerance.
    while (dL > toler) {
        z.normalise();
        q = z;
        z = A * q;
        // Update the eigenvalue estimate and its relative change.
        real Lo = eigenValue;
        eigenValue = dot(q, z);
        dL = fabs(Lo - eigenValue) / eigenValue;
        iter++;
        if (iter >= maxIter) {
            break;
        }
    }
    eigenVector = q;
    return iter;
}
} // end namespace SpaMtrix
