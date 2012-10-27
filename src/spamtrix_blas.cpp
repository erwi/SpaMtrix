#include <setup.h>
#include <spamtrix_blas.h>

namespace SpaMtrix
{

real errorNorm2(const IRCMatrix& A,
		const Vector& x,
		const Vector& b)
{
  /*!
   * Calculates error nor for system Ax = b, using 
   * error = (Ax-b).(Ax-b)
   */
  real error(-1);
  
  Vector temp(x.getLength() );
  multiply(A, x, temp);  	// temp = Ax
  temp-= b;
  
  error = dot(temp,temp); 
  return error;
}

real errorNorm2(const TDMatrix& A,
		const Vector& x,
		const Vector& b)
{
  /*!
   * Calculates error nor for system Ax = b, using 
   * error = (Ax-b).(Ax-b)
   */
  real error(-1);
  
  Vector temp(x.getLength() );
  multiply(A, x, temp);  	// temp = Ax
  temp-= b;
  
  error = dot(temp,temp); 
  return error;
}

} // end namespace SpaMtrix