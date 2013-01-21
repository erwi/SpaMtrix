#ifndef MATRIXMAKER_H
#define MATRIXMAKER_H

#include <vector>
#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <cstring>
#include <cmath>

#include <spamtrix_setup.hpp>
#include <spamtrix_ircmatrix.hpp>
#include <spamtrix_fleximatrix.hpp>

namespace SpaMtrix
{

class MatrixMaker
{
    /*!
   * Class for specifying matrix sparsity pattern and creating the corresponding matrix.
   */
    idx nRows;
    idx nCols;

    //std::vector< std::list< idx > > nonZeros;
    FlexiMatrix nz;
    MatrixMaker(){}
public:

    MatrixMaker(const idx nRows, const idx nCols);

    idx calcNumNonZeros() const;
    void addNonZero(const idx row, const idx col, const real val = 0.0);
    void expandBlocks(const idx numExp = 1);

    void poisson5Point(); // CREATES A 5 POINT POISSON FINITE DIFFERENCES TEST MATRIX


    IRCMatrix getIRCMatrix();
    void makeSparseMatrix(IRCMatrix &A); //CONVERTS A SPARSITY PATTERN TO A SPARSE MATRIX
};
} // end namespace SpaMtrix

#endif

