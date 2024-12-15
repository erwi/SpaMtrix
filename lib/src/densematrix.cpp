#include <stdlib.h>
#include <iostream>
#include <stdio.h>

#include <spamtrix_densematrix.hpp>
#include <spamtrix_exception.hpp>

using namespace SpaMtrix;
DenseMatrix::DenseMatrix(idx numRows, idx numCols):
        numRows(numRows),
        numCols(numCols) {
  values = std::vector<real>(numRows * numCols, 0.0);
  if (values.size() !=  numRows * numCols) {
    throw SpaMtrixException("error in " + std::string(ERROR_LOCATION) + " could not allocate dense matrix");
  }
}
DenseMatrix::DenseMatrix(const DenseMatrix &other):
        numRows(other.numRows),
        numCols(other.numCols) {
  values = other.values;
}

DenseMatrix::~DenseMatrix() = default;

void DenseMatrix::setAllValuesTo(const real v) {
  for (auto &value : values) {
    value = v;
  }
}

real &DenseMatrix::operator()(const idx row, const idx col) {
  /*!
    Row/column access operator
    */
#ifdef DEBUG
  assert(row < numRows);
    assert(col < numCols_);
#endif
  const idx ind = col * numRows + row;
  return values[ind];
}

real DenseMatrix::operator()(const idx row, const idx col) const {
  /*!
    Row/column access operator
    */
#ifdef DEBUG
  assert(row < numRows);
    assert(col < numCols_);
#endif
  const idx ind = col * numRows + row;
  return values[ind];
}


void DenseMatrix::print() const {
  for (idx r = 0 ; r < numRows ; r++) {
    for (idx c = 0  ; c < numCols ; c++) {
      printf("%1.2f ", values[c * numRows + r]);
    }
    printf("\n"); fflush(stdout);
  }
}
