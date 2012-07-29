#ifndef SMOOTHERS_H
#define SMOOTHERS_H
#include <setup.h>
#include <ircmatrix.h>
#include <vector.h>

/*!
 * THE JACOBI, GAUSS-SEIDEL AND SUCCESSIVE-OVER-RELAXATION (SOR) SOLVERS/SMOOTHERS
 * ARE DECLARED AND DEFINED HERE
 */
void jacobi(const IRCMatrix& A,
	    Vector &x,
	    const Vector &b,
	    const idx &maxIter)
{
  /*!
   * Performs maxIter Jacobi iterations to solve Ax=b
   */
  Vector x_old( x.getLength() );
  for (idx k = 0 ; k < maxIter ; ++k) // FOR maxIter
  {
    x_old = x;
    for (idx i = 0 ; i < A.getNumRows() ; i++ ) // FOR ROWS
    {
      real rho(0.0);
      real diagonal(1.0);
      // PERFORM MATRIX-VECTOR MULTIPLICATION, NEGLECTING THE DIAGONAL TERM
      const idx row_start = A.rows[i];
      const idx row_end = A.rows[i+1];
      for (idx j = row_start ; j < row_end ; ++j)
      {
    const idx col = A.cvPairs[j].ind;
	if (col != i ) // IF NOT DIAGONAL ENTRY
	{
	  rho += A.cvPairs[j].val * x_old[j];
	}//
	else
	{
	  diagonal = A.cvPairs[j].val; // STORE DIAGONAL VALUE FOR LATER USE
	}
      }// end for j
      x[i] = (b[i] - rho) / diagonal; 
    }// end for i
    
  }// end for k
  
  
}

void gauss_seidel(const IRCMatrix& A, 
		  Vector &x, 
		  const Vector &b, 
		  idx& maxIter)
{
  /*!
   * makes maxIter Gauss-Seidel iterations to solve Ax=b
   */

  for (idx k = 0 ; k < maxIter ; ++k) // FOR ITERATIONS
  {
    for (idx i = 0 ; i < A.getNumRows() ; ++i)
    {
      real rho = 0;
      for ( idx j = 0 ; j < A.getNumRows() ; ++j)
      {
	rho+= A.sparse_get(i,j)*x[j]*(i!=j);	// ONLY NON-DIAGONALS
      }// end for j
      x[i] = 1/( A.sparse_get(i,i) ) * (b[i] - rho );
    }// end for i
  }//end for k
}



//===============================
// SUCCESSIVE OVER RELAXATION
void SOR(const IRCMatrix& A, 
  Vector &x, 
  const Vector &b, 
  idx maxIter)
{
  
  real omega = 1.5;
  
  for (idx k = 0 ; k < maxIter ; ++k)
  {
    for (idx i = 0 ; i < A.getNumRows() ; ++i)
    {
      real rho(0);
      real diagonal(1.0);
      idx col_start = A.rows[i];
      idx col_end = A.rows[i+1];
      
      for (idx j = col_start ; j < col_end ; ++j)
      {
    if ( i ==A.cvPairs[j].ind ) // IF DIAGONAL
	{
	  diagonal = omega / A.cvPairs[j].val;
	}
	else // OFF-DIAGONAL
	{
      rho += A.cvPairs[j].val*x[A.cvPairs[j].ind]; // rho = rho + A[i,j]*x[j]
	}
      }// end for j
      x[i] = (1.0-omega)*x[i] + diagonal*(b[i] - rho );
      
    }// end for i
  }//end for k
}



#endif


