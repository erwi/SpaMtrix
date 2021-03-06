#include <spamtrix_matrixmaker.hpp>
#include <spamtrix_ircmatrix.hpp>
#include <algorithm>
#include <cassert>
#include <cmath>

namespace SpaMtrix {

    MatrixMaker::MatrixMaker(const idx nRows, const idx nCols) :
            nRows(nRows),
            nCols(nCols),
            nz(nRows, nCols) {
    }

    MatrixMaker::~MatrixMaker() = default;

    /*!
      Adds storage for a non-zero at location row, col. If value of non-zero is
      known, it is added too, otherwise zero is assumed.
    */
    void MatrixMaker::addNonZero(const idx row, const idx col, const real val) {
        nz.addNonZero(row, col, val);
    }

    /*!Returns total number of nonzeros allocated so far.*/
    idx MatrixMaker::calcNumNonZeros() const {
        return nz.calcNumNonZeros();
    }

    /*!
    * Makes a 5 point finite differences (2D) poisson test matrix.
    * The grid spacing is assumed to be unity, resulting in a main
    * diagonal with value 4, and all off-diagonals with values -1.
    *
    * The FD grid is assumed to have an equal number of rows and columns n,
    * where n = sqrt( side length ) of the built matrix A.
    *
    * Example:
    *          To build matrix that corresponds to a 10 x 10 FD grid,
    *          create a MatrixMake object:
    *             1. MatrixMaker mm(100, 100).
    *             2. mm.poisson5Point().
    *             3. IRCMatrix A = mm.getIRCMatrix();
    */
    void MatrixMaker::poisson5Point() {

        // MAKE SURE MATRIX IS SQUARE AND NON-ZERO SIZED
        assert(nRows == nCols);
        assert(nRows);
        idx n = sqrt(nRows); // FINITE DIFFERENCES GRID SIDE LENGTH
        // SET MAIN DIAGONALS TO 4
        for (idx i = 0; i < nRows; i++)
            addNonZero(i, i, 4);
        // SET OFF-DIAGONALS
        for (idx i = 0; i < nRows; i++) {
            idx row = i / n; // ROW OF i'th NODE
            idx col = i % n; // COLUMN OF i'th NODE
            // RIGHT DERIVATIVE
            if (col < n - 1)
                addNonZero(i, i + 1, -1);
            // LEFT DERIVATIVE
            if (col > 0)
                addNonZero(i, i - 1, -1);
            // UP DERIVATIVE
            if (row > 0)
                addNonZero(i - n, i, -1);
            if (row < n - 1)
                addNonZero(i + n, i, -1);
        }// end off-diagonals
    }// end void poisson5Point

    void MatrixMaker::expandBlocks(const idx numExp) {
        /*!
          Expands existing sparsity pattern by a factor of numExp+1. E.g., if numExp is 2, the resulting matrix size
          will be 3 times the original: [a] -> |aaa|
                                               |aaa|
                                               |aaa|
        */
        if (!numExp)
            return;
        // EXPAND TO RIGHT: a -> aa...a
        // FOR EACH ROW
        for (idx r = 0; r < nRows; ++r) {
            idx numC = nz.nonZeros[r].size();
            for (idx e = 1; e <= numExp; ++e) { // for number of expansions
                for (idx c = 0; c < numC; ++c) {
                    IndVal iv = nz.nonZeros[r][c];
                    iv.ind += e * nCols;
                    nz.nonZeros[r].push_back(iv);
                    //nz.addNonZero(r,iv); // <-- SHOULD USE THIS INSTEAD
                    nz.numDim2 = nz.numDim2 > (iv.ind + 1) ? nz.numDim2 : (iv.ind + 1);
                }// end for c
            }//end for e
        }// end for r
        nCols *= (numExp + 1);
        // EXPAND DOWN: aaa -> aaa
        //                     aaa
        //                     aaa
        for (idx e = 1; e <= numExp; ++e) {
            for (idx r = 0; r < nRows; ++r) {
                std::vector<IndVal>::iterator rb = nz.nonZeros[r].begin();
                std::vector<IndVal>::iterator re = nz.nonZeros[r].end();
                nz.nonZeros.push_back(std::vector<IndVal>(rb, re));
            }// end for r
        }//end for e
        nRows *= (numExp + 1);
        nz.numDim1 = nz.nonZeros.size();
    }// end void expandBlocks

    IRCMatrix MatrixMaker::getIRCMatrix() {
        /*!
        * Converts sparsity pattern to IRCMatrix \n
        * Use makeSparseMatrix instead!
        */
        return std::move(IRCMatrix(nz));
    }// end getIRCMatrix

    void MatrixMaker::makeSparseMatrix(IRCMatrix &A) {
        /*!
         * converts a sparsity pattern to a Sparse matrix
         */
        A.copyFrom(nz);
    }// end makeSparseMatrix
} // end namespace SpaMtrix
