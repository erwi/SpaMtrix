// Copyright (c) 2026 Eero
// SPDX-License-Identifier: MIT
#ifndef FLEXIMATRIX_H
#define FLEXIMATRIX_H
#include <vector>
#include <spamtrix_setup.hpp>
#include <spamtrix_ircmatrix.hpp>

namespace SpaMtrix {
  class IRCMatrix;

  /**
   * @brief Flexible sparse matrix storage using rows of index-value pairs.
   */
  class FlexiMatrix {
      size_t numCols_ = 0;
      std::vector<std::vector<IndVal> > nonZeros;
    public:

    /** @brief Construct an empty matrix. */
    FlexiMatrix(): numCols_(0) {}

    /**
     * @brief Construct from a sparse IRC matrix.
     *
     * @param A Source matrix.
     */
    FlexiMatrix(const IRCMatrix &A);
    /** @brief Destroy the matrix. */
    virtual ~FlexiMatrix();

    /**
     * @brief Add storage for a non-zero entry.
     *
     * @param row Row index.
     * @param col Column index.
     * @param value Initial value.
     */
    void addNonZero(size_t row, size_t col, real value = 0.0);

    /**
     * @brief Count the number of stored non-zero entries.
     *
     * @return Number of stored entries.
     */
    [[nodiscard]] idx calcNumNonZeros() const;

    /**
     * @brief Return the value at a matrix position.
     *
     * @param row Row index.
     * @param col Column index.
     * @return Stored value, or zero if the entry is not present.
     */
    [[nodiscard]] real getValue(size_t row, size_t col) const;
    /** @brief Get the number of rows. */
    [[nodiscard]] size_t getNumRows() const { return nonZeros.size(); };
    /** @brief Get the number of columns. */
    [[nodiscard]] size_t getNumCols() const { return numCols_; };

    /**
     * @brief Access a row of stored non-zero entries.
     *
     * @param row Row index.
     * @return Mutable row storage.
     */
    [[nodiscard]] std::vector<IndVal>& row(size_t row);
    /**
     * @brief Access a row of stored non-zero entries.
     *
     * @param row Row index.
     * @return Const row storage.
     */
    [[nodiscard]] const std::vector<IndVal>& row(size_t row) const;

    /**
     * @brief Set a matrix entry, inserting storage if needed.
     *
     * @param dim1 Row index.
     * @param dim2 Column index.
     * @param val Value to assign.
     */
    void setValue(const idx dim1, const idx dim2, const real val);
    /** @brief Print the matrix to stdout. */
    void print() const;

    /**
     * @brief Check whether a position is stored as non-zero.
     *
     * @param dim1 Row index.
     * @param dim2 Column index.
     * @param val Optional pointer to receive the stored value.
     * @return True if the entry exists.
     */
    [[nodiscard]] bool isNonZero(const idx dim1, const idx dim2, real *&val);
  };
} // end namespace SpaMtrix










#endif // FLEXIMATRIX_H

