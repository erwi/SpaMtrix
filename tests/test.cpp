#define CATCH_CONFIG_MAIN
#include <catch.h>
#include <spamtrix_matrixmaker.hpp>
using namespace SpaMtrix;

TEST_CASE("MatrixMaker creates IRCMatrix", "[MatrixMaker]") {
    REQUIRE(true);
    // Create sparse matrix like:
    //    ######
    //    .#####
    //    ..####
    //    ...###
    //    ....##

    idx numRows = 5;
    idx numCols = 6;
    int numAdded = 0;
    MatrixMaker matrixMaker(numRows, numCols);
    // Create sparse matrix with non-zeros in upper diagonal only
    for (idx r = 0; r < numRows; r++) {
        for (idx c = 0; c < numCols; c++) {
            if (r <= c) {
                matrixMaker.addNonZero(r, c);
                numAdded ++;
                REQUIRE(numAdded == matrixMaker.calcNumNonZeros());
            }
        }
    }

    IRCMatrix matrix = matrixMaker.getIRCMatrix();
    REQUIRE(matrix.getNumRows() == numRows);
    REQUIRE(matrix.getNumCols() == numCols);
    REQUIRE(20 == matrix.getnnz());

    for (idx r = 0; r < numRows; r++) {
        for (idx c = 0; c < numCols; c++) {
            if (r <= c) {
                REQUIRE(matrix.isNonZero(r, c));
            } else {
                REQUIRE(!matrix.isNonZero(r, c));
            }
        }
    }
}

TEST_CASE("Expand block along matrix diagonal") {
    // Define sparsity pattern for 3x3 matrix like:
    // a.d
    // .b.
    // ..c

    const double a = 1, b = 2, c = 3, d = 13;

    MatrixMaker mm(3, 3);
    mm.addNonZero(0, 0, a);
    mm.addNonZero(1, 1, b);
    mm.addNonZero(2, 2, c);
    mm.addNonZero(0, 2, d);

    // Expand the sparsity pattern along diagonal to the 9x9 matrix
    // a.d......
    // .b.......
    // ..c......
    // ...a.d...
    // ....b....
    // .....c...
    // ......a.d
    // .......b.
    // ........c

    mm.expandDiagonal(3);

    auto M = mm.getIRCMatrix();

    REQUIRE(M.getNumRows() == 9);
    REQUIRE(M.getNumCols() == 9);
    REQUIRE(M.getnnz() == 12);

    for (int r = 0; r < 9; r++) {
        for (int c = 0; c < 9; c++) {
            if (r == c) {
                // diagonal positions
                REQUIRE(M.getValue(r, c) == r % 3 + 1);
            }
            else if ((r == 0 && c == 2) || (r == 3 && c == 5) || (r == 6 && c == 8)) {
                // off-diagonal non-zero
                REQUIRE(M.getValue(r, c) == d);
            } else {
                // everything else should be zero sparse
                REQUIRE(M.isNonZero(r, c) == false);
            }
        }
    }
}