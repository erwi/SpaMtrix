#ifndef DENSEMATRIX_H
#define DENSEMATRIX_H

#include <vector>
#include <spamtrix_setup.hpp>
namespace SpaMtrix{
class DenseMatrix
{
    /*!
     * DenseMatrix is a column ordered dense matrix
     */
    std::vector<real> values;
    idx numRows;
    idx numCols;

    DenseMatrix();
    void setAllValuesTo(const real v = 0);
public:
    DenseMatrix(const idx numRows, const idx numCols);
    DenseMatrix(const DenseMatrix &other);
    ~DenseMatrix();
    real& operator()(const idx row, const idx col);
    real operator()(const idx row, const idx col) const;
    idx getNumRows()const {return numRows;}
    idx getNumCols()const {return numCols;}
    
    void print() const;
};
} // end namespace SpaMtrix
#endif // DENSEMATRIX_H
