// Copyright (c) 2026 Eero
// SPDX-License-Identifier: MIT
#include <iostream>
#include <spamtrix_matrixmaker.hpp>
#include <spamtrix_ircmatrix.hpp>
#include <spamtrix_vector.hpp>
#include <spamtrix_cholesky.hpp>
#include <spamtrix_blas.hpp>

using namespace SpaMtrix;
int main()
{
  std::cout <<"Solving linear system Ax=b using Cholesky decomposition"<<std::endl;
  // Create the sparse matrix.
  MatrixMaker mm(5,5);
  
  // Row 1.
  mm.addNonZero(0,0, 24);  mm.addNonZero(0,2,6);
  // Row 2.
  mm.addNonZero(1,1,8);  mm.addNonZero(1,2,2);
  // Row 3.
  mm.addNonZero(2,0,6);  mm.addNonZero(2,1,2);
  mm.addNonZero(2,2,8);  mm.addNonZero(2,3,-6);
  mm.addNonZero(2,4,2);
  // Row 4.
  mm.addNonZero(3,2,-6);  mm.addNonZero(3,3,24);
  // Row 5.
  mm.addNonZero(4,2,2);  mm.addNonZero(4,4,8);

  IRCMatrix A = mm.getIRCMatrix();

  std::cout<<"Matrix A : "<<std::endl;
  A.print();

  Cholesky M(A);
  // Create vectors x and b for the system Ax=b.
  Vector x(A.getNumCols());
  Vector b(x);
  b[0] = 1.0;   // Initial condition.

  // Solve the system using forward and back substitution.
  std::cout << "Solving Ax=b...";
  M.solve(x,b);
  std::cout << "OK" << std::endl;
  // Print the numerical error magnitude.
  real e = sqrt(errorNorm2(A,x,b));
  std::cout<< "error is " << e << std::endl;

  return 0;
}
