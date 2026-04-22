// Copyright (c) 2026 Eero
// SPDX-License-Identifier: MIT
#ifndef POWERMETHOD_H
#define POWERMETHOD_H
#include <spamtrix_blas.hpp>
#include <spamtrix_setup.hpp>
namespace SpaMtrix{
/**
 * @brief Estimate the dominant eigenpair using power iteration.
 *
 * @param A Input matrix.
 * @param eigenValue Output dominant eigenvalue estimate.
 * @param eigenVector Input initial guess; output eigenvector estimate.
 * @param toler Requested relative tolerance.
 * @param maxIter Maximum number of iterations.
 * @return Number of iterations performed.
 */
idx powerMethod(const SpaMtrix::IRCMatrix &A,
                 real &eigenValue,
                 SpaMtrix::Vector &eigenVector,
                 real &toler,
                 unsigned int maxIter = MAX_INDEX
                );

} // end namespace SpaMtrix
#endif
