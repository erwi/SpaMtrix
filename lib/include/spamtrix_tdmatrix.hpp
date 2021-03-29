#ifndef TDMATRIX_H
#define TDMATRIX_H

#include <iostream>
#include <stdlib.h>
#include <string.h>

#include <spamtrix_setup.hpp>
#include <spamtrix_vector.hpp>
namespace SpaMtrix
{
/*!
* Declares a class for tridiagonal sparse matrix
*/

class TDMatrix{

  real* upper;
  real* diagonal;
  real* lower;
  idx size;

  bool isValidIndex(const idx row, const idx col) const;

public:
  TDMatrix(const idx size);
  virtual ~TDMatrix();
  idx getnnz() const {return 3*size-2;}
  idx getNumRows() const {return size;}
  idx getNumCols() const {return size;}
  void sparse_set(const idx row, const idx col, const real val);
  void sparse_add(const idx row, const idx col, const real val);
  real sparse_get(const idx row, const idx col)const;

  void solveAxb(Vector& x, const Vector& b) const;
  void print(const char* name = NULL) const;

  //=============================================
  // FRIEND FUNCTIONS
  //
  friend void multiply(const TDMatrix& A, const Vector& x, Vector& b); // defined in spamtrix_blas.h
};

} // end namespace SpaMtrix
#endif

