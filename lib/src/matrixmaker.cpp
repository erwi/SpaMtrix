#include <spamtrix_matrixmaker.hpp>
#include <spamtrix_ircmatrix.hpp>
#include <cassert>
#include <cmath>

namespace SpaMtrix {

    MatrixMaker::MatrixMaker(const idx nRows, const idx nCols) :
            nRows(nRows),
            nCols(nCols), nz()  { }

    MatrixMaker::~MatrixMaker() = default;

    /*!
      Adds storage for a non-zero at location row, col. If value of non-zero is
      known, it is added too, otherwise zero is assumed.
    */
    void MatrixMaker::addNonZero(const idx row, const idx col, const real val) {
      assert(row < nRows);
      assert(col < nCols);
      nz.addNonZero(row, col, val);
    }

    /*!Returns total number of nonzeros allocated so far.*/
    idx MatrixMaker::calcNumNonZeros() const {
        return nz.calcNumNonZeros();
    }

    void MatrixMaker::poisson5Point() {

        // MAKE SURE MATRIX IS SQUARE AND NON-ZERO SIZED
        assert(nRows == nCols);
        assert(nRows);
        idx n = sqrt(nRows); // FD grid side length
        assert(n * n == nRows); // verify that a valid grid length is provided

        // Set main diagonals to 4
        for (idx i = 0; i < nRows; i++) {
          addNonZero(i, i, 4);
        }
        // Set off-diagonal terms to -1
        for (idx i = 0; i < nRows; i++) {
            idx row = i / n; // ROW OF i'th NODE
            idx col = i % n; // COLUMN OF i'th NODE
            // RIGHT DERIVATIVE
            if (col < n - 1) {
              addNonZero(i, i + 1, -1);
            }
            // LEFT DERIVATIVE
            if (col > 0) {
              addNonZero(i, i - 1, -1);
            }
            // UP DERIVATIVE
            if (row > 0) {
              addNonZero(i - n, i, -1);
            }
            if (row < n - 1) {
              addNonZero(i + n, i, -1);
            }
        }
    }

    void MatrixMaker::expandBlocks(const idx numExp) {

      if (!numExp)
        return;
      // Expand each row to right: a -> aa...a
      for (idx r = 0; r < nRows; ++r) {
        auto &rowNonZeros = nz.row(r);
        const idx numC = rowNonZeros.size(); //nonZeros[r].size();
        for (idx e = 1; e <= numExp; ++e) { // for number of expansions
          for (idx c = 0; c < numC; ++c) {
            real val = rowNonZeros[c].val;
            idx col = rowNonZeros[c].ind;

            nz.addNonZero(r, col + e * nCols, val);
            //IndVal iv = nz.nonZeros[r][c];
            //iv.ind += e * nCols;
            //nz.nonZeros[r].push_back(iv);
            //nz.addNonZero(r,iv); // <-- SHOULD USE THIS INSTEAD
            //nz.numDim2 = nz.numDim2 > (iv.ind + 1) ? nz.numDim2 : (iv.ind + 1);
          }// end for c
        }//end for e
      }// end for r
      nCols *= (numExp + 1);
      // Expand rows down: aaa -> aaa
      //                     aaa
      //                     aaa
      for (idx e = 1; e <= numExp; ++e) {
        for (idx r = 0; r < nRows; ++r) {
          for (auto &nonZero : nz.row(r)) {
            idx col = nonZero.ind;
            real val = nonZero.val;
            nz.addNonZero(r + e * nRows, col, val);
          }
        }
      }
      nRows *= (numExp + 1);
    }

    IRCMatrix MatrixMaker::getIRCMatrix() {
        /*!
        * Converts sparsity pattern to IRCMatrix \n
        * Use makeSparseMatrix instead!
        */
        return std::move(IRCMatrix(nz));
    }// end getIRCMatrix

    IRCMatrix* MatrixMaker::newIRCMatrix() {
      return new IRCMatrix(nz);
    }

    void MatrixMaker::makeSparseMatrix(IRCMatrix &A) {
        /*!
         * converts a sparsity pattern to a Sparse matrix
         */
        A.copyFrom(nz);
    }

    void MatrixMaker::expandDiagonal(idx numExp) {
        if (numExp <= 1) {
            return;
        }
        const size_t numRowsInitial = nz.getNumRows();
        const size_t numColsInitial = nz.getNumCols();
        for (idx i = 1; i < numExp; i++) {
            for (size_t row = 0; row < numRowsInitial; row++) {
                const size_t newRowIdx = (i * numRowsInitial) + row;
                for (const IndVal &nonZero : nz.row(row)) {
                    const size_t newColIdx = (i * numColsInitial) + nonZero.ind;
                    nz.addNonZero(newRowIdx, newColIdx, nonZero.val);
                }
            }
        }
        nRows = nz.getNumRows();
        nCols = nz.getNumCols();
    }
}