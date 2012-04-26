#ifndef MATRIXMAKER_H
#define MATRIXMAKER_H
#include <setup.h>
#include <ircmatrix.h>
#include <vector>
#include <list>
#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <cstring>
using std::vector;
using std::list;
using std::cout;
using std::endl;
class MatrixMaker
{
  /*!
   * Class for specifying matrix sparsity pattern and creating the corresponding matrix.
   */
  idx nRows;
  idx nCols;
  
  std::vector< std::list< idx > > nonZeros; 
  
public:
  MatrixMaker();
  //~MatrixMaker(); // IMPLICIT IS OK
  void setMatrixSize(const idx nRows, const idx nCols);
  void addNonZero(const idx row, const idx col); 

  void clear();
  
  IRCMatrix getIRCMatrix();
  
};


#endif

