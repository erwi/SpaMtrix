#include <tdmatrix.h>

TDMatrix::TDMatrix(const idx size):
upper(NULL),
diagonal(NULL),
lower(NULL),
size(size)
{
  upper = new real [size-1];
  diagonal = new real [size];
  lower = new real [size-1];
  
  if ( (!upper) || (!diagonal) || (!lower) )
  {
    cout <<"error in " << __func__ << "could not allocate for matrix of size " << size << endl;
    exit(1);
  }
 
}

TDMatrix::~TDMatrix()
{
  size = 0;
  delete [] upper; 
  upper = NULL;
  delete [] lower;
  lower = NULL;
  delete [] diagonal;
  diagonal = NULL;
}

bool TDMatrix::isValidIndex(const idx row, const idx col) const
{
  /*!
   * Checks for out of bounds indexing. TODO
   */

  
  return true;
}


void TDMatrix::sparse_set(const idx row, const idx col, const real val)
{
  #ifdef DEBUG
    assert(isValidIndex(row,col);
  #endif
      
     if (row==col)
     {
      diagonal[row] = val;
      return;
     }
    else if (col > row) 
    {
      upper[col-1] = val;
      return;
    }
    else
      lower[row-1] = val;
    
}

void TDMatrix::sparse_add(const idx row, const idx col, const real val) 
{
  #ifdef DEBUG
    assert(isValidIndex(row,col);
  #endif
    
     
     if (row==col)
     {
      diagonal[row] += val;
      return;
     }
    else if (col > row) 
    {
      upper[col-1] += val;
      return;
    }
    else
      lower[row-1] += val;
    
}

real TDMatrix::sparse_get(const idx row, const idx col) const
{
  #ifdef DEBUG
    assert(isValidIndex(row,col);
  #endif
    
    if (row==col)
      return diagonal[row];
    else if (col > row) 
      return upper[col-1];
    else
      return lower[row-1];
      
}

void TDMatrix::solveAxb(Vector& x, const Vector& b) const
{
  /*!
   * Solves Ax=b using Tridiagonal Matrix Algorithm (TDMA) a.k.a. Thomas algorithm.
   * See wikipedia http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
   */
  
  /*
   
     a - sub-diagonal (means it is the diagonal below the main diagonal) -- indexed from 1..n-1
    * d - the main diagonal
    * c - sup-diagonal (means it is the diagonal above the main diagonal) -- indexed from 0..n-2
    * v - right part
    * x - the answer
    */
  
  real *a = this->lower;
  real *c = this->upper;
  
  real *d = new real[size]; // ALLOCATE TEMPORARY DIAGONAL
  real *v = new real[size]; // TEMPORARY RHS VECTOR
  
// COPY TEMPORARY VARIABLES
  for (idx i = 0 ; i < size ; i++)
  {
    d[i] = diagonal[i];
    v[i] = b[i];
  }
  

  for ( idx i = 1 ; i < size; ++i)
  {
    if ( d[i-1] == 0.0 )
    {
      cout << "error in " << __func__ << " diagonal "<<i<< "is zero"<< endl;
      exit(1);
    }
    double m = a[i] / d[i-1];
    
    d[i] -= m*c[i-1];
    v[i] -= m*v[i-1];
  }
  
  x[size-1] = v[size-1] / d[size-1];
  
  for (idx i = size - 1; i > 0 ; --i)
    x[i-1] = (v[i-1]- c[i-1] * x[i] ) / d[i-1];
  
  delete [] d;
  delete [] v;
}

void TDMatrix::print()
{
  /*
  cout << "Tridiagonal matrix size: "<< size << endl;
  for (idx i = 0 ; i < size - 1 ; i++)
    cout << "l,d,u " << i <<" = " << lower[i] << "," << diagonal[i] << "," << upper[i] << endl;
  
  cout << "  d   " << size-1<< " = " << diagonal[size-1] << endl;
  */
  
  for (idx r = 0 ; r < size ; r++)
  {
    for idx c = 0 ; c < size ; c++)
    {
      
    }
  }
    
}


