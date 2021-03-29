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