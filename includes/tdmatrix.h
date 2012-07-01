#ifndef TDMATRIX_H
#define TDMATRIX_H
#include <setup.h>
#include <iostream>
#include <stdlib.h>
#include <vector.h>
#include <string.h>
/*!
* Declares a class for tridiagonal sparse matrix
*/
using std::cout;
using std::endl;
class TDMatrix{
  real* upper;
  real* diagonal;
  real* lower;
  idx size;
  
  bool isValidIndex(const idx row, const idx col) const;
  
public:
  TDMatrix(const idx size);
  ~TDMatrix();
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


#endif

