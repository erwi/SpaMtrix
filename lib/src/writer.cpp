#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <spamtrix_writer.hpp>
#include <spamtrix_vector.hpp>
#include <spamtrix_ircmatrix.hpp>
namespace SpaMtrix {

Writer::~Writer() {
}

bool Writer::writeCSV(const std::string &filename,
                      const Vector &data,
                      idx numC) {
    int numD = data.getLength();
    if (numC == 0) {
        numC = numD;
    }
    // Open the file for writing.
    std::fstream file;
    file.open(filename.c_str() , std::fstream::out);
    if (!file.is_open()) {
        return false;
    }
    // Write all data.
    idx cc = 0;
    for (int i = 0; i < numD; i++) {
        file << data[i];
        cc++;
        if (cc < numC) { // IF NOT END OF ROW ...
            file << ",";
        } else { // ROW END REACHED - NEW LINE
            file << "\n";
            cc = 0;
        }
    }
    file.close();
    return true;
} // end writeCSV


bool Writer::writeMatrixMarket(const std::string &filename,
                               const IRCMatrix &A) {
    std::fstream file;
    // Open the file for writing.
    file.open(filename.c_str(), std::fstream::out);
    if (!file.is_open()) {
        return false;
    }
        // Write the header.
    file << "% MatrixMarket matrix coordinate real general" << std::endl;
    file << "% NOTE: Matrix market indexing is 1-based!" << std::endl;
        // Write the matrix size descriptor.
    file << A.getNumRows() << " " << A.getNumCols() << " " << A.getnnz()
         << std::endl;
        // Write matrix data.
    for (idx r = 0; r < A.getNumRows(); r++) {
        for (idx c = 0; c < A.getNumCols(); c++) {
            real val;
            if (A.isNonZero(r, c, val)) {
                file << r + 1 << "\t" << c + 1 << "\t" << val << std::endl;
            }
        }
    }
    return true;
}// end bool writeMatrixMarket
} // end namespace SpaMtrix
