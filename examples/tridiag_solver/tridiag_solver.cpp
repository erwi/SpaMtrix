// Copyright (c) 2026 Eero
// SPDX-License-Identifier: MIT
#include <iostream>
#include <math.h>

// SpaMtrix includes.
#include <spamtrix_ircmatrix.hpp>
#include <spamtrix_vector.hpp>
#include <spamtrix_matrixmaker.hpp>
#include <spamtrix_blas.hpp>
#include <spamtrix_diagpreconditioner.hpp>
#include <spamtrix_tdmatrix.hpp>

// Solve a 1D Poisson finite-difference problem.
using std::cout;
using std::endl;
using namespace SpaMtrix;
int main(int nargs, char* args[])
{
  // Construct the 1D finite-difference matrix.
  unsigned int np = 10;
  
  if (nargs > 1){
    np = atoi( args[1] );
  }
  
  
  cout <<"Solving : Ax=b"<<endl;
  cout << "creating tridiagonal matrix of size:" << np <<"x" <<np<<"...";
  

  // Fill in the matrix values.
  //		    | 2  -1     |
  // A = 1/(h^2) *  | -1  2  -1 |
  //		    |    -1   2 |		  
  //
  // Using h = 1.
  
  TDMatrix tdm(np); // Create the tridiagonal matrix.

  for (unsigned int i = 0 ; i < np ; i++ ){
    tdm.sparse_set(i,i,2.0);  // Diagonal.
    if ( i > 0 ){
      tdm.sparse_set(i,i-1, -1.0); // Sub-diagonal.
    }
    if ( i < np - 1){
      tdm.sparse_set(i, i+1 , -1); // Super-diagonal.
    }
  }
  cout<<"OK"<<endl;  
  // Create the unknown vector x with fixed values 1 and -1 at both ends.
  Vector x(np);
  x[0] = 1.0;
  x[np-1] = -1.0;
   
  // Right-hand-side vector b.
  Vector b(np);
  // Apply boundary conditions.
  multiply(tdm,x,b); 	// b = Ax;
 
  scale(-1.0, b); // WANT TO SOLVE Ax = -b, SO MULTIPLY BY -1 HERE
  
  // Modify matrix rows and columns for known nodes.
  // First node.
  tdm.sparse_set(0, 0, 1.0);
  tdm.sparse_set(0, 1, 0.0);
  tdm.sparse_set(1, 0, 0.0);
  // Last node.
  tdm.sparse_set(np-1, np-1, 1.0);
  tdm.sparse_set(np-2, np-1, 0.0);
  tdm.sparse_set(np-1, np-2, 0.0);
  
 // Print the vectors if the problem is small.
  if (np <= 10){
    cout << "R.H.S. vector b:"<<endl;
    b.print("b");
    tdm.print("A");
  }
    
  // SOLVE Ax=b
  tdm.solveAxb(x,b);
   
  cout << "error after TDM solver : "<< sqrt(errorNorm2(tdm,x,b)) << endl;
  
  // Restore fixed boundary node values.
  x[0] = 1.0;
  x[np-1] = -1.0;
  
  // Print the result if the problem is small.
  if (np <= 10){
    cout << "solution vector x:"<< endl;
    x.print("x");
  }
 
  return 0;
}
