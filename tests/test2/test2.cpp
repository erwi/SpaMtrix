#include <ircmatrix.h>
#include <vector.h>
#include <matrixmaker.h>
#include <spamtrix_blas.h>

#include <diagpreconditioner.h>
#include <jacobipreconditioner.h>
#include <sorpreconditioner.h>
#include <pcg.h>
#include <cg.h>
#include <iostream>
#include <math.h>
#include <tdmatrix.h>
// TEST 2. SOLVES 1D POISSON FINITE DIFFERENCES PROBLEM
using std::cout;
using std::endl;
int main(int nargs, char* args[])
{
  // CONSTRUCT 1D FINITE DIFFERENCES MATRIX
  unsigned int np = 10;
  
  if (nargs > 1)
    np = atoi( args[1] );
  
  
  cout <<"Solving : Ax=b"<<endl;
  cout << "creating tridiagonal matrix of size:" << np <<"x" <<np<<"...";
  

  // FILL IN MATRIX VALUES
  //		    | 2  -1     |
  // A = 1/(h^2) *  | -1  2  -1 |
  //		    |    -1   2 |		  
  //
  // USING h = 1:
  
  TDMatrix tdm(np); // MAKE TRIDIAGONAL MATRIX

  for (unsigned int i = 0 ; i < np ; i++ )
  {
    tdm.sparse_set(i,i,2.0);  // DIAGONAL
    if ( i > 0 )
    {
      tdm.sparse_set(i,i-1, -1.0); // SUB-DIAGONAL
    }
    if ( i < np - 1)
    {
      tdm.sparse_set(i, i+1 , -1); // SUP-DIAGONAL
    }
  }
  cout<<"OK"<<endl;  
  // CREATE UNKNOWN VECTOR x, WITH FIXED VALUES 1, -1 AT BOTH ENDS
  Vector x(np);
  x[0] = 1.0;
  x[np-1] = -1.0;
   
  // R.H.S VECTOR b
  Vector b(np);
  // ===================================== 
  // APPLY BOUNDARY CONDITIONS
  // =====================================
  multiply(tdm,x,b); 	// b = Ax;
 
  scale(-1.0, b); // WANT TO SOLVE Ax = -b, SO MULTIPLY BY -1 HERE
  
  // MODIFY MATRIX COLUMNS/ROWS FOR KNOWN NODES
  // FIRST NODE
  tdm.sparse_set(0, 0, 1.0);
  tdm.sparse_set(0, 1, 0.0);
  tdm.sparse_set(1, 0, 0.0);
  // LAST NODE
  tdm.sparse_set(np-1, np-1, 1.0);
  tdm.sparse_set(np-2, np-1, 0.0);
  tdm.sparse_set(np-1, np-2, 0.0);
  
 // PRINT IF SMALL PROBLEM
  if (np <= 10)
  {
    cout << "R.H.S. vector b:"<<endl;
    b.print("b");
    tdm.print("A");
  }
    
  // SOLVE Ax=b
  tdm.solveAxb(x,b);
   
  cout << "error after TDM solver : "<< sqrt(errorNorm2(tdm,x,b)) << endl;
  
  // RESTORE FIXED BOUNDARY NODE VALUES
  x[0] = 1.0;
  x[np-1] = -1.0;
  
  // PRINT RESULT IF SMALL PROBLEM
  if (np <= 10)
  {
    cout << "solution vector x:"<< endl;
    x.print("x");
  }
 
  return 0;
}