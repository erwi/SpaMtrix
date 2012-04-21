#include <setup.h>
#include <spamtrix_blas.h>

void multiply(const IRCMatrix& A,
	     const Vector& x,
	     Vector& b)
{
  /*!
   * MATRIX VECTOR MULTIPLICATION Ax=b
   */
  
#ifdef DEBUG
  assert(A.getNumCols() == x.getLength() );
#endif
  // IF b AND x ARE NOT REFERENCES TO SAME VECTOR
  if (&b != &x)
  {
/*
    (if (b.getLength() != x.getLength() ) // RESIZE b IF NECESSARY
    {
      b.resize( x.getLength() )
    }
    
  */    
    // FOR EACH ROW
    for (idx i = 0 ; i < A.getNumRows() ; i++)
    {
      // FOR EACH COLUMN
      real r(0);
      for (idx j = 0 ; j < A.getNumCols() ; j++)
      {
	r += A.sparse_get(i,j)*x[j];
      }//j
      b[i] = r;
      
    }//
    
    
  }
  
  else // b OVERWRITES X
  {
    // WILL NEED TO USE TEMPORARY VECTOR VARIABLE HERE
    std::cout << " option not supported yet in " << __func__<<std::endl;
    exit(1);
  }
  
  
}

// VECTOR SCALAR PRODUCT v2 = a*v1
void multiply(const Vector& v1, const real a, Vector& v2)
{
#ifdef DEBUG
assert(v1.getLength() == v2.getLength() );
#endif
  for (idx i = 0 ; i < v1.getLength() ; ++i)
    v2[i] = v1[i]*a;
}

real dot(const Vector& v1, const Vector& v2)
{
  /*!
   * Calculates dot product between two vectors v1 and v2
   */
#ifdef DEBUG
  assert(v1.getLength() == v2.getLength() );
#endif
  real d(0.0);
  for (idx i = 0 ; i < v1.getLength() ; ++i)
    d+= v1[i] * v2[i];
  
  return d;
}

