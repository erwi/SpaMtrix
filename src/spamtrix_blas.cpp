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

    // FOR EACH ROW
#ifdef USES_OPENMP
#pragma omp parallel for
#endif
    for (idx i = 0 ; i < A.getNumRows() ; i++)
    {
      // FOR EACH COLUMN
      real r(0);
      idx row_start = A.rows[i];
      idx row_end   = A.rows[i+1];
      for (idx j = row_start ; j < row_end ; j++)
      {
	idx col = A.cvPairs[j].col;
	r += A.cvPairs[j].val * x[col];
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
#pragma omp parallel for reduction(+:d)
  for (idx i = 0 ; i < v1.getLength() ; ++i)
    d+= v1[i] * v2[i];
  
  return d;
}

