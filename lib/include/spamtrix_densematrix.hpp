#ifndef DENSEMATRIX_H
#define DENSEMATRIX_H
#include <vector>
#include <spamtrix_setup.hpp>

namespace SpaMtrix {
  /*!
  * DenseMatrix is a column ordered dense matrix
  */
  class DenseMatrix {
    std::vector<real> values;
    idx numRows;
    idx numCols;

    DenseMatrix();
  public:
    DenseMatrix(idx numRows, idx numCols);
    DenseMatrix(const DenseMatrix &other);
    virtual ~DenseMatrix();

    void setAllValuesTo(real v);
    real& operator()(idx row, idx col);
    real operator()(idx row, idx col) const;
    [[nodiscard]] idx getNumRows() const { return numRows; }
    [[nodiscard]] idx getNumCols() const { return numCols; }

    void print() const;
  };
}
#endif // DENSEMATRIX_H
