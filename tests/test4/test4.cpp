/*!
 * Test file for cholesky decomposition
 * Using a symmetric positive definite matrix A:
 * 24 0  6  0  0
 *  0 8  2  0  0
 *  6 2  8 -6  2
 *  0 0 -6 24  0
 *  0 0  2  0  8
 */


#include <iostream>
#include <matrixmaker.h>
#include <ircmatrix.h>
#include <cholesky.h>
#include <spamtrix_blas.h>
int main()
{
  std::cout <<"Solving linear system Ax=b using Cholesky decomposition"<<std::endl;
// CREATE EMPTY SPARSE MATRIX
  MatrixMaker mm;
  mm.setMatrixSize(5,5);
  // ROW 1
  mm.addNonZero(0,0);  mm.addNonZero(0,2);
  // ROW 2
  mm.addNonZero(1,1);  mm.addNonZero(1,2);
  // ROW 3
  mm.addNonZero(2,0);  mm.addNonZero(2,1);
  mm.addNonZero(2,2);  mm.addNonZero(2,3);
  mm.addNonZero(2,4);
  // ROW 4
  mm.addNonZero(3,2);  mm.addNonZero(3,3);
  // ROW 5
  mm.addNonZero(4,2);  mm.addNonZero(4,4);

  IRCMatrix A = mm.getIRCMatrix();

// SET VALUES IN SPARSE MATRIX A
  A.sparse_set(0,0, 24); A.sparse_set(0,2, 6);
  A.sparse_set(1,1,8); A.sparse_set(1,2,2);
  A.sparse_set(2,0,6); A.sparse_set(2,1,2);
  A.sparse_set(2,2,8); A.sparse_set(2,3,-6);
  A.sparse_set(2,4,2);
  A.sparse_set(3,2,-6); A.sparse_set(3,3,24);
  A.sparse_set(4,2,2); A.sparse_set(4,4,8);
  
  std::cout<<"Matrix A : "<<std::endl;
  A.print();

  Cholesky M(A);
  //M.print();

  // CREATE VECTORS x AND b FOR SYSTEM Ax=b
  Vector x(A.getNumCols());
  Vector b(x);
  b[0] = 1.0;   // INITIAL CONDITION

  // SOLVE SYSTEM USING FORWARD/BACK-SUBSTITUTION
  std::cout << "Solving Ax=b...";
  M.solve(x,b);
  std::cout << "OK" << std::endl;
  // PRINT NUMERICAL ERROR MAGNITUDE
  real e = sqrt(errorNorm2(A,x,b));
  std::cout<< "error is " << e << std::endl;

  return 0;
}
