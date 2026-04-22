#ifndef MATRIXMAKER_H
#define MATRIXMAKER_H
#include <vector>
#include <spamtrix_setup.hpp>
#include <spamtrix_fleximatrix.hpp>

namespace SpaMtrix {
class IRCMatrix;

/**
 * @brief Helper for building sparse matrix sparsity patterns.
 */
class MatrixMaker {
    idx nRows;          // NUMBER OF ROWS
    idx nCols;          // NUMBER OF COLUMSN
    FlexiMatrix nz;     // TEMPORARY "FLEXIBLE" SPARSE MATRIX DATASTUCTURE
    MatrixMaker(){}
public:
    /**
     * @brief Construct a matrix maker for the requested size.
     *
     * @param nRows Number of rows.
     * @param nCols Number of columns.
     */
    MatrixMaker(const idx nRows, const idx nCols);
    /** @brief Destroy the helper. */
    virtual ~MatrixMaker();
    /** @brief Count the currently stored non-zero entries. */
    idx calcNumNonZeros() const;
    /**
     * @brief Add a non-zero position to the sparsity pattern.
     *
     * @param row Row index.
     * @param col Column index.
     * @param val Initial value.
     */
    void addNonZero(const idx row, const idx col, const real val = 0.0);

  /**
   * @brief Expand each stored entry into a square block pattern.
   *
   * If @p numExp is 2, each entry becomes a $3\times 3$ block.
   *
   * @param numExp Expansion factor.
   */
  void expandBlocks(const idx numExp = 1);
    /**
     * @brief Expand the diagonal pattern to a larger diagonal block matrix.
     *
     * @param numExp Expansion factor.
     */
    void expandDiagonal(idx numExp);

    /**
     * @brief Build a 5-point finite-difference Poisson test matrix.
     *
     * The grid spacing is assumed to be unity.
     */
    void poisson5Point();

    /**
     * @brief Create an identity matrix sparsity pattern.
     *
     * The matrix must be square.
     */
    void identity();

    /**
     * @brief Materialize the current sparsity pattern as an IRC matrix.
     *
     * @return The constructed sparse matrix.
     */
    IRCMatrix getIRCMatrix();
    /**
     * @brief Allocate a new IRC matrix on the heap.
     *
     * @return Newly allocated matrix.
     */
    IRCMatrix* newIRCMatrix();
    /**
     * @brief Populate an existing IRC matrix from the current pattern.
     *
     * @param A Matrix to populate.
     */
    void makeSparseMatrix(IRCMatrix &A);
};
}

#endif

