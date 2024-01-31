#ifndef MATRIXMAKER_H
#define MATRIXMAKER_H
#include <vector>
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
    /**
     * Expands non-zeroes along matrix diagonal to size numExp.
     * For example the N-by-N block on non-zeroes # expanded with numExp=3 becomes
     * a 3N-by-3N <br>
     *     #.. <br>
     *     .#. <br>
     *     ..# <br>
     *
     */
    void expandDiagonal(idx numExp);
    void poisson5Point(); // CREATES A 5 POINT POISSON FINITE DIFFERENCES TEST MATRIX
    IRCMatrix getIRCMatrix();
    /** Creates a new new IRCMatrix on the heap and returns pointer to it.
     * Caller is responsible of deleting the pointer
     */
    IRCMatrix* newIRCMatrix();
    void makeSparseMatrix(IRCMatrix &A); //CONVERTS A SPARSITY PATTERN TO A SPARSE MATRIX
};// end class MatrixMaker
} // end namespace SpaMtrix

#endif

