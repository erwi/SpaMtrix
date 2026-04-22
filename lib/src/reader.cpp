// Copyright (c) 2026 Eero
// SPDX-License-Identifier: MIT
#include <iostream>
#include <ios>
#include <sstream>
#include <string>

#include <spamtrix_reader.hpp>
#include <spamtrix_ircmatrix.hpp>
#include <spamtrix_matrixmaker.hpp>
#include <spamtrix_exception.hpp>

namespace SpaMtrix {
const char *Reader::FILE_OPEN_ERROR_STRING = "Could not open file : ";
const char *Reader::FILE_FORMAT_ERROR_STRING = "Bad format reading file : ";
IRCMatrix Reader::readMatrixMarket(const std::string &filename) {
    // Open the file.
    std::ifstream file;
    file.open(filename);
    if (!file.is_open()) {
        throw SpaMtrixException(FILE_OPEN_ERROR_STRING + filename + " at " + std::string(ERROR_LOCATION));
    }
    // Skip comment lines at the top of the file.
    std::string line = "%";
    while (line.at(0) == '%') {
        std::getline(file, line);
    }
    // The next line contains the row, column, and non-zero counts.
    idx numRows = 0, numCols = 0, numNZ = 0;
    if (!(std::stringstream(line) >> numRows >> numCols >> numNZ)) {
        throw (SpaMtrixException(FILE_FORMAT_ERROR_STRING + filename + " at " + std::string(ERROR_LOCATION)));
    }
    MatrixMaker mm(numRows, numCols);
    // Read all non-zero entries and add them to the matrix maker.
    for (idx i = 0 ; i < numNZ ; ++i) {
        std::getline(file, line);
        idx row, col;
        real val;
        if (!(std::stringstream(line) >> row >> col >> val)) {
            throw (SpaMtrixException(FILE_FORMAT_ERROR_STRING + filename + " at " + std::string(ERROR_LOCATION)));
        }
        mm.addNonZero(row - 1, col - 1, val); // TAKE INTO ACCOUNT 1-BASED INDEXING
    }
    file.close();
    return mm.getIRCMatrix();
} // end readMatrixMarket
}// end namespace SpaMtrix
