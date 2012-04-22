#include <pcg.h>

bool pcg(const IRCMatrix& A,
	 const Vector& b,
	 Vector& x,
	 const DiagPreconditioner& M,
	 idx& maxIter,
	 real& toler)
{
/*!
 * Preconditioned conjugate gradient solver for solving Ax=b.
 * Returns true if converged to solution with tolerance toler
 * in maxIter iterations. Total number of iterations needed is returned
 * in maxIter and tolerance achieved in toler 
 */
  Vector temp1( b.getLength() );
  Vector temp2( b.getLength() );
  
  // INITIAL RESIDUAL r = b - Ax
  multiply (A,x,temp1);
  Vector r = b;
  r-=temp1;
  
  // d = M^(-1) r
  Vector d(r);
  M.applyToVector(d);
  
  real delta_n = dot(r,d);
  real delta_0 = delta_n;
  unsigned int i = 0;
  
  while ( (i < maxIter) && (delta_n > toler*toler*delta_0) )
  {
    
  }
  
  return false;
}