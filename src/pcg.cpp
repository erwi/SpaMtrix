#include <pcg.h>

bool pcg(const IRCMatrix& A,
	 const Vector& b,
	 Vector& x,
	 const Preconditioner& M,
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
  Vector d(b.getLength() );
  M.solveMxb(d,r);
  
  real delta_n = dot(r,d);
  real delta_0 = delta_n;
  unsigned int i = 0;
  
  while ( (i < maxIter) && (delta_n > toler*toler*delta_0) )
  {
    // q = Ad
    multiply(A,d,temp1); // temp1 = q
    
    // CALCULATE alpha
    real alpha = delta_n / dot(d,temp1);
 
    // x = x + alpha*d
    multiply(d,alpha,temp2); // temp2 = alpha*d
    x+= temp2;
    
    if ( i % 50 == 0 )
    {
      // r = b - Ax
      r = b;
      multiply(A,x, temp2 );
      r-= temp2;
    }
    else
    {
      // r = r - alpha*q
      multiply(temp1,alpha, temp2); // temp2 = alpha*q
      r-= temp2;
    }

    // s = M^(-1)*r
    //temp1 = r;
    //M.applyToVector(temp1); // temp1 = s
    M.solveMxb(temp1,r);
    
    
    real delta_old = delta_n;
    delta_n = dot(r,temp1);
    
    real beta = delta_n / delta_old;
    
    // d = s + beta*d
    multiply(d,beta,d); // d*beta
    d += temp1; 	// + s
    
    i++;
  }
    
  bool converged( i < maxIter);
  maxIter = i;
  toler = delta_n;
  return converged;
}