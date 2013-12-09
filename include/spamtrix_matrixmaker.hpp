#ifndef MATRIXMAKER_H
#define MATRIXMAKER_H

#include <vector>
#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <cstring>
#include <cmath>

#include <spamtrix_setup.hpp>
#include <spamtrix_fleximatrix.hpp>

namespace SpaMtrix{
class IRCMatrix;

class MatrixMaker{
/*!
* Class for specifying matrix sparsity pattern and creating the corresponding sparse matrix data structure.
*/
    idx nRows;          // NUMBER OF ROWS
    idx nCols;          // NUMBER OF COLUMSN
    FlexiMatrix nz;     // TEMPORARY "FLEXIBLE" SPARSE MATRIX DATASTUCTURE
    MatrixMaker(){}
public:
    MatrixMaker(const idx nRows, const idx nCols);
    virtual ~MatrixMaker();
    idx calcNumNonZeros() const;
    void addNonZero(const idx row, const idx col, const real val = 0.0);
    void expandBlocks(const idx numExp = 1);
    void poisson5Point(); // CREATES A 5 POINT POISSON FINITE DIFFERENCES TEST MATRIX
    IRCMatrix getIRCMatrix();
    void makeSparseMatrix(IRCMatrix &A); //CONVERTS A SPARSITY PATTERN TO A SPARSE MATRIX
};// end class MatrixMaker
} // end namespace SpaMtrix

#endif

