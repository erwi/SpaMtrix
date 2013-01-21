#ifndef DENSEMATRIX_H
#define DENSEMATRIX_H

#include <vector>

#include <spamtrix_setup.hpp>
class DenseMatrix
{
    //real *values;
    std::vector<real> values;
    idx numRows;
    idx numCols;

    DenseMatrix();
    void setAllValuesTo(const real v = 0);
public:
    DenseMatrix(const idx numRows, const idx numCols);
    ~DenseMatrix();
    real& operator()(const idx row, const idx col);
    real operator()(const idx row, const idx col) const;
    void print() const;
};

#endif // DENSEMATRIX_H
