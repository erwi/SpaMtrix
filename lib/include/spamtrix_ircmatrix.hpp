// Copyright (c) 2026 Eero
// SPDX-License-Identifier: MIT
#ifndef IRCMATRIX
#define IRCMATRIX

#include <spamtrix_setup.hpp>
namespace SpaMtrix
{
// forward decarations
class FlexiMatrix;
class Vector;

/**
 * @brief Sparse matrix stored in Interleaved Row Compressed form.
 */
class IRCMatrix{
protected:
    idx* rows;      // ROW COUNTER
    IndVal* cvPairs;// COLUMN-VALUE PAIRS
    idx nnz;        // NUMBER OF NON-ZEROS
    idx numRows, numCols;

    IndVal & find(const idx row, const idx col);
    [[nodiscard]] const IndVal& find(const idx row, const idx col) const;
public:

  /** @brief Construct an empty matrix. */
    IRCMatrix();
  /**
   * @brief Construct a matrix from raw sparse storage.
   *
   * @param numRows Number of rows.
   * @param numCols Number of columns.
   * @param nnz Number of non-zero entries.
   * @param rows Row pointer array.
   * @param cvPairs Column-value pairs.
   */
    IRCMatrix(  const idx numRows, const idx numCols,
                const idx nnz,
                idx * const rows, IndVal *const cvPairs);
  /** @brief Copy-construct a matrix. */
    IRCMatrix(const IRCMatrix &m);
  /** @brief Move-construct a matrix. */
    IRCMatrix(IRCMatrix &&m);
  /**
   * @brief Construct from a flexible sparse matrix.
   *
   * @param M Source matrix.
   */
  IRCMatrix(const FlexiMatrix &M);
  /** @brief Copy-assign a matrix. */
    IRCMatrix& operator=(const IRCMatrix& m);
  /** @brief Move-assign a matrix. */
    IRCMatrix& operator=(IRCMatrix&& m);
  /** @brief Fill all stored entries with a scalar value. */
    IRCMatrix& operator=(const real &s);
  /** @brief Assign from a flexible sparse matrix. */
    IRCMatrix& operator=(const FlexiMatrix &m);
  /** @brief Destroy the matrix. */
    virtual ~IRCMatrix();
    //================================================
  /** @brief Release all owned storage. */
    void clear();
  /**
   * @brief Copy data from a flexible sparse matrix.
   *
   * @param A Source matrix.
   */
    void copyFrom(const FlexiMatrix& A);
  /** @brief Get the number of non-zero entries. */
  idx getnnz()const;
  /** @brief Get the number of rows. */
  idx getNumRows()const;
  /** @brief Get the number of columns. */
  idx getNumCols() const;

  /**
   * @brief Set a sparse entry.
   *
   * @param row Row index.
   * @param col Column index.
   * @param val Value to store.
   */
    void sparse_set(const idx row, const idx col , const real val );
  /**
   * @brief Add a value to a sparse entry.
   *
   * @param row Row index.
   * @param col Column index.
   * @param val Value to add.
   */
    void sparse_add(const idx row, const idx col , const real val );
  /**
   * @brief Read a sparse entry.
   *
   * @param row Row index.
   * @param col Column index.
   * @return Stored value, or zero if the entry is not present.
   */
    real sparse_get(const idx row, const idx col ) const;

  /**
   * @brief Read a matrix entry.
   *
   * @param row Row index.
   * @param col Column index.
   * @return Matrix value, or zero if the entry is not stored.
   */
    real getValue(const idx row, const idx col) const;

  /**
   * @brief Get a pointer to a stored value.
   *
   * @param row Row index.
   * @param col Column index.
   * @return Pointer to the stored value, or nullptr if absent.
   */
    [[nodiscard]] real* getValuePtr(const idx row, const idx col);

  /**
   * @brief Get a const pointer to a stored value.
   *
   * @param row Row index.
   * @param col Column index.
   * @return Pointer to the stored value, or nullptr if absent.
   */
    [[nodiscard]] real* getValuePtr(const idx row, const idx col) const;

  /**
   * @brief Multiply the matrix by a vector.
   *
   * @param x Input vector.
   * @return The product $A x$.
   */
  Vector operator*(const Vector& x) const;
  /**
   * @brief Scale the matrix in place.
   *
   * @param s Scaling factor.
   */
  void operator*=(const real &s);
  /**
   * @brief Return a scaled copy of the matrix.
   *
   * @param s Scaling factor.
   * @return The scaled matrix.
   */
  const IRCMatrix operator*(const real &s) const;

    /**
   * @brief Add another matrix, optionally scaled.
   *
   * The sparsity patterns must match.
   *
   * @param other Matrix to add.
   * @param scalar Scaling factor applied to @p other.
     */
    void add(const IRCMatrix& other, const real& scalar = 1.0);

  /**
   * @brief Check whether storage exists at a position.
   *
   * @param row Row index.
   * @param col Column index.
   * @return True if the entry is stored.
   */
    [[nodiscard]] bool isNonZero(const idx row, const idx col) const;
  /**
   * @brief Check whether storage exists and copy the value.
   *
   * @param row Row index.
   * @param col Column index.
   * @param val Output value.
   * @return True if the entry is stored.
   */
    [[nodiscard]] bool isNonZero(const idx row, const idx col, real& val) const;
    //===============================================
    // FRIEND FUNCTIONS THAT REQUIRE ACCESS TO PRIVATE DATA FOR PERFORMANCE.
    friend void multiply(const IRCMatrix& A, const Vector& x, Vector& b);
    friend real multiply_dot(const IRCMatrix& A, const Vector& x, Vector& b);
    friend class FlexiMatrix;
    //================================================
    // DEBUG FUNCTIONS
    /** @brief Print the sparsity pattern for debugging. */
    void spy()const;
    /** @brief Print the matrix to stdout. */
    void print() const;
};

inline const IRCMatrix operator*(const real& s, const IRCMatrix &M)
{
    return M*s;
}

// IMPLEMENTATIONS OF INLINDED METHODS - NO NEW DECLARATIONS BELOW THIS
inline idx IRCMatrix::getNumCols() const { return numCols; }
inline idx IRCMatrix::getNumRows() const { return numRows; }
inline idx IRCMatrix::getnnz() const { return nnz;}


} // end namespace SpaMtrix
#endif
