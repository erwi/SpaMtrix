#include <cg.h>
#include <spamtrix_blas.h>
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
  Vector temp1( b.getLength() );
  Vector temp2( b.getLength() );
  multiply(A,x,temp1);		// temp = Ax
  r-= temp1;
  
  // d = r
  Vector d(r);
  
  real delta_n = dot(r,r);
  real delta_0 = delta_n;
  
  // START MAIN LOOP
  idx i = 0;
  
  while ( (i < maxIter) && (delta_n > toler*toler*delta_0 ) )
  {
  
      multiply(A,d,temp1); // temp = Ad (= q)
      real a = delta_n / dot(d,temp1);
      
      // x = x + a*d
      multiply(d,a,temp2);
      x+=temp2;
      
      if (i % 50 == 0 )
      {
	
	// r = b-Ax
	r = b;
	multiply(A,x, temp2);
	r-=temp2;
      }
      else
      {
	// r = r-a*q
	multiply(temp1, a, temp2);
	r-=temp2;
      }
            
      real delta_old = delta_n;
      delta_n =  dot(r,r);
      real Beta = delta_n / delta_old;
      
      // d = r + Beta*d
      multiply(d,Beta,temp2); 	// Beta*r
      d = r;
      d+= temp2;
      
      i++;
  
  }
  
  // NUMBER OF ITERATIONS PERFORMED
  bool converged(i < maxIter);
  maxIter = i; 
  
  toler = delta_n;
  
  
  return converged;
  
}
