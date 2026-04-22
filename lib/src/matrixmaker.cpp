// Copyright (c) 2026 Eero
// SPDX-License-Identifier: MIT
#include <spamtrix_matrixmaker.hpp>
#include <spamtrix_ircmatrix.hpp>
#include <cassert>
#include <cmath>

namespace SpaMtrix {

    MatrixMaker::MatrixMaker(const idx nRows, const idx nCols) :
            nRows(nRows),
            nCols(nCols), nz()  { }

    MatrixMaker::~MatrixMaker() = default;

    void MatrixMaker::addNonZero(const idx row, const idx col, const real val) {
      assert(row < nRows);
      assert(col < nCols);
      nz.addNonZero(row, col, val);
    }

    idx MatrixMaker::calcNumNonZeros() const {
        return nz.calcNumNonZeros();
    }

    void MatrixMaker::poisson5Point() {

        // Ensure the matrix is square and non-empty.
        assert(nRows == nCols);
        assert(nRows);
        idx n = sqrt(nRows); // FD grid side length
        assert(n * n == nRows); // verify that a valid grid length is provided

        // Set the main diagonal to 4.
        for (idx i = 0; i < nRows; i++) {
          addNonZero(i, i, 4);
        }
        // Set the off-diagonal terms to -1.
        for (idx i = 0; i < nRows; i++) {
            idx row = i / n; // ROW OF i'th NODE
            idx col = i % n; // COLUMN OF i'th NODE
            // Right neighbour.
            if (col < n - 1) {
              addNonZero(i, i + 1, -1);
            }
            // Left neighbour.
            if (col > 0) {
              addNonZero(i, i - 1, -1);
            }
            // Upper neighbour.
            if (row > 0) {
              addNonZero(i - n, i, -1);
            }
            if (row < n - 1) {
              addNonZero(i + n, i, -1);
            }
        }
    }

    void MatrixMaker::identity() {
        assert(nRows == nCols);
        for (idx i = 0; i < nRows; i++) {
            addNonZero(i, i, 1);
        }
    }

    void MatrixMaker::expandBlocks(const idx numExp) {

      if (!numExp)
        return;
      // Expand each row to the right.
      for (idx r = 0; r < nRows; ++r) {
        auto &rowNonZeros = nz.row(r);
        const idx numC = rowNonZeros.size();
        for (idx e = 1; e <= numExp; ++e) {
          for (idx c = 0; c < numC; ++c) {
            real val = rowNonZeros[c].val;
            idx col = rowNonZeros[c].ind;

            nz.addNonZero(r, col + e * nCols, val);
          }
        }
      }
      nCols *= (numExp + 1);
      // Expand rows downward.
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
        return std::move(IRCMatrix(nz));
    }// end getIRCMatrix

    IRCMatrix* MatrixMaker::newIRCMatrix() {
      return new IRCMatrix(nz);
    }

    void MatrixMaker::makeSparseMatrix(IRCMatrix &A) {
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