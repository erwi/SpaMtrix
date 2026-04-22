// Copyright (c) 2026 Eero
// SPDX-License-Identifier: MIT
#ifndef TDMATRIX_H
#define TDMATRIX_H

#include <iostream>
#include <stdlib.h>
#include <string.h>

#include <spamtrix_setup.hpp>
#include <spamtrix_vector.hpp>
namespace SpaMtrix
{
/**
 * @brief Tridiagonal sparse matrix.
 */
class TDMatrix{

  real* upper;
  real* diagonal;
  real* lower;
  idx size;

  bool isValidIndex(const idx row, const idx col) const;

public:
  /**
   * @brief Construct a tridiagonal matrix of the requested size.
   *
   * @param size Matrix dimension.
   */
  TDMatrix(const idx size);
  /** @brief Destroy the matrix. */
  virtual ~TDMatrix();
  /** @brief Get the number of non-zero entries. */
  idx getnnz() const {return 3*size-2;}
  /** @brief Get the number of rows. */
  idx getNumRows() const {return size;}
  /** @brief Get the number of columns. */
  idx getNumCols() const {return size;}
  /**
   * @brief Set a sparse entry.
   *
   * @param row Row index.
   * @param col Column index.
   * @param val Value to store.
   */
  void sparse_set(const idx row, const idx col, const real val);
  /**
   * @brief Add a value to a sparse entry.
   *
   * @param row Row index.
   * @param col Column index.
   * @param val Value to add.
   */
  void sparse_add(const idx row, const idx col, const real val);
  /**
   * @brief Read a sparse entry.
   *
   * @param row Row index.
   * @param col Column index.
   * @return Stored value, or zero if the entry is not present.
   */
  real sparse_get(const idx row, const idx col)const;

  /**
   * @brief Solve $A x = b$.
   *
   * @param x Solution vector to populate.
   * @param b Right-hand-side vector.
   */
  void solveAxb(Vector& x, const Vector& b) const;
  /**
   * @brief Print the matrix to stdout.
   *
   * @param name Optional label.
   */
  void print(const char* name = NULL) const;

  //=============================================
  // FRIEND FUNCTIONS
  friend void multiply(const TDMatrix& A, const Vector& x, Vector& b); // defined in spamtrix_blas.h
};

} // end namespace SpaMtrix
#endif

