#ifndef DENSEMATRIX_H
#define DENSEMATRIX_H
#include <vector>
#include <spamtrix_setup.hpp>

namespace SpaMtrix {
  /**
   * @brief Column-major dense matrix.
   */
  class DenseMatrix {
    std::vector<real> values;
    idx numRows;
    idx numCols;

    DenseMatrix();
  public:
    /**
     * @brief Construct a dense matrix with the given dimensions.
     *
     * @param numRows Number of rows.
     * @param numCols Number of columns.
     */
    DenseMatrix(idx numRows, idx numCols);
    /**
     * @brief Copy-construct a dense matrix.
     *
     * @param other Source matrix.
     */
    DenseMatrix(const DenseMatrix &other);
    /** @brief Destroy the matrix. */
    virtual ~DenseMatrix();

    /**
     * @brief Set every stored value to the same scalar.
     *
     * @param v Value to write.
     */
    void setAllValuesTo(real v);
    /**
     * @brief Access a mutable matrix entry.
     *
     * @param row Row index.
     * @param col Column index.
     * @return Reference to the selected entry.
     */
    real& operator()(idx row, idx col);
    /**
     * @brief Access a matrix entry.
     *
     * @param row Row index.
     * @param col Column index.
     * @return Value at the selected entry.
     */
    real operator()(idx row, idx col) const;
    /** @brief Get the number of rows. */
    [[nodiscard]] idx getNumRows() const { return numRows; }
    /** @brief Get the number of columns. */
    [[nodiscard]] idx getNumCols() const { return numCols; }

    /** @brief Print the matrix to stdout. */
    void print() const;
  };
}
#endif // DENSEMATRIX_H
