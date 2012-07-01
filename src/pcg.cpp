#include <pcg.h>
#include <math.h>
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
 * in maxIter and tolerance achieved in toler.  
 */
  
// INITIAL RESIDUAL
  Vector r = b;
  Vector Av(x.getLength() );
  multiply(A,x,Av);
  r-= Av;

  // z = M^(-1)r
  Vector z(x.getLength() );
  M.solveMxb(z,r);
  
  // p = z
  Vector p(z);
  idx k(0);
  bool converged = false;
  real error_new = dot(r,z);
   
  while ( ( k < maxIter ) &&
	  (error_new > toler*toler) )
  {
    
    // Av = A*p
    multiply(A,p,Av);
    const real alpha = error_new / dot (Av, p);
    

    // UPDATE SOLUTION x = x + alpha*p

    axpy(alpha, p, x);
    
    // UPDATE RESIDUAL r = r - alpha*Ap

    axpy(-1*alpha, Av, r);
    
       
    M.solveMxb(z,r);
    const real error_old = error_new;
    error_new = dot(r,z);
    const real beta = error_new / error_old ;
    // p = beta*p + z
    aypx(beta,p,z);
        
    k++;
  }
  converged = sqrt(error_new) <= toler;
  
  maxIter = k;
  toler = sqrt(error_new);
  
  return converged;
}