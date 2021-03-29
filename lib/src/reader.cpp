#include <iostream>
#include <ios>
#include <sstream>
#include <string>

#include <spamtrix_reader.hpp>
#include <spamtrix_ircmatrix.hpp>
#include <spamtrix_matrixmaker.hpp>

namespace SpaMtrix {
const char *Reader::FILE_OPEN_ERROR_STRING = "Could not open file : ";
const char *Reader::FILE_FORMAT_ERROR_STRING = "Bad format reading file : ";
IRCMatrix Reader::readMatrixMarket(const std::string &filename) {
    /*!
       Reads sparse matrix from a matrix market formated (text) file.\n
       This assumes that no lines may start with spaces/tabs/other whitespace \n
       characters. Also, no empty lines are allowed.\n \n
       Exceptions are thrown if file could not be opened or if the file format \n
       is incorrect.
     */
    // ATTEMPT TO OPEN FILE
    std::ifstream file;
    file.open(filename);
    if (!file.is_open()) {
        throw std::ios_base::failure(FILE_OPEN_ERROR_STRING + filename);
    }
    // FILE MAY START WITH LINES OF COMMENTS.
    // COMMENTED LINES BEGIN WITH THE '%'-CHARACTER
    std::string line = "%";
    while (line.at(0) == '%') {
        std::getline(file, line);
    }
    // THE VARIABLE 'line' NOW CONTAINS NUMBER OF ROWS, COLUMNS AND NONZEROS
    idx numRows = 0, numCols = 0, numNZ = 0;
    if (!(std::stringstream(line) >> numRows >> numCols >> numNZ)) {
        throw (std::ios_base::failure(FILE_FORMAT_ERROR_STRING + filename));
    }
    MatrixMaker mm(numRows, numCols);
    // READ ALL NONZERO ENTRIES AND ADD TO MATRIXMAKER
    for (idx i = 0 ; i < numNZ ; ++i) {
        std::getline(file, line);
        idx row, col;
        real val;
        if (!(std::stringstream(line) >> row >> col >> val)) {
            throw (std::ios_base::failure(FILE_FORMAT_ERROR_STRING + filename));
        }
        mm.addNonZero(row - 1, col - 1, val); // TAKE INTO ACCOUNT 1-BASED INDEXING
    }
    file.close();
    return mm.getIRCMatrix();
} // end readMatrixMarket
}// end namespace SpaMtrix
