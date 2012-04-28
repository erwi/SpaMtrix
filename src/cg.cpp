#include <cg.h>
//#include <spamtrix_blas.h>
bool cg(const IRCMatrix& A,
	const Vector& b,
	Vector& x,
	idx& maxIter,
	real& toler)
{
// CONJUGATE GRADIENT SOLVER
// SOLVES Ax = b
// RETURNS TRUE IF CONVERGED TO SOLUTION IN maxIter ITERATIONS
// TOTAL NUMBER OF ITERATIONS IS RETURNED IN maxIter
// TOLERANCE ACHIEVED IS RETURNED IN toler
  
  
  // r = b - Ax
  Vector r = b;
  Vector Av(r.getLength() ); // TEMPORARY VECTOR FOR MATRIX*VEC 
  multiply(A, x, Av);
  r-=Av;
  // p = r;
  Vector p = r;
  
  // START MAIN LOOP
  idx i = 0;
  
  real error_old = dot(r,r);
  cout << "initial error :" << error_old << endl;
  while ( (i < maxIter) && ( error_old > toler*toler ) )
    //(delta_n > //toler*toler*delta_0 ) )
  {
      // CALCULATE alpha
      multiply(A,p,Av); 
      real alpha = error_old / dot(p,Av);
      
      // x = x + alpha*d
      axpy(alpha,p,x);
        
      //r = r-alpha*Ad
      axpy( -1.0*alpha , Av, r);
            
      const real error_new = dot(r,r);
      if ( error_new < toler*toler )
      {
	error_old = error_new;
	break;
      }
      
      printf("iter %i , error %e\n", (int) i, error_new);
      const real error_ratio = error_new / error_old;
      // p = error_ratio*p + r
      aypx(error_ratio, p, r);
      error_old = error_new;
      i++;
  }
  
  // NUMBER OF ITERATIONS PERFORMED
  bool converged(i < maxIter);
  maxIter = i; 
  
  toler = error_old;
  
  
  return converged;
  
}
