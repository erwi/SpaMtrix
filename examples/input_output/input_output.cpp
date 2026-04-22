// Copyright (c) 2026 Eero
// SPDX-License-Identifier: MIT
#include <iostream>
#include <string>

#include <spamtrix_reader.hpp>
#include <spamtrix_writer.hpp>
#include <spamtrix_matrixmaker.hpp>
#include <spamtrix_ircmatrix.hpp>
int main(int nargs, char* args[] ){
    idx size = 9;
    if (nargs>1)
        size = atoi(args[1]);
        
    // 1. Create the test matrix.
    SpaMtrix::MatrixMaker mm(size,size);
    mm.poisson5Point();
    SpaMtrix::IRCMatrix A = mm.getIRCMatrix();
    A.spy();
    
    // 2. Write the test matrix to a file.
    std::string filename("testoutput.txt");
    SpaMtrix::Writer::writeMatrixMarket(filename,A);
    A.clear();
    
    // 3. Read the test matrix from the file.
    SpaMtrix::IRCMatrix B =
    SpaMtrix::Reader::readMatrixMarket(filename);
    B.print();
    
    return 0;
}
