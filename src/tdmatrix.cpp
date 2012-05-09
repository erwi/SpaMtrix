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

  if ( (row >= size ) || (col>=size) )
    return false;
  
  if (row == col ) return true; 	// diagonal
  if (row == col+1) return true;	// lower
  if (col == row+1) return true;	// upper
  
  return false;	// otherwise not tridiagonal
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

void TDMatrix::solveAxb(Vector& x, Vector& b) //const
{
  /*!
   * Solves Ax=b using Tridiagonal Matrix Algorithm (TDMA) a.k.a. Thomas algorithm.
   * See wikipedia http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
   * WARNING: Currently this modified the matrix!!
   */
  idx N = this->size;	// LENGTH OF VECTOR
  
  // FORWARD LOOP - ELIMINATE SUB-DIAGONAL
  // MODIFY FIRST-ROW COEFFS
  upper[0] = upper[0] / diagonal[0];
  b[0] = b[0] / diagonal[0];
  for (idx i = 1 ; i < N-1 ; ++i)
  {
    real temp = diagonal[i] - lower[i-1]*upper[i-1];
    upper[i] = upper[i] / temp;
    b[i] = (b[i] - lower[i-1]*b[i-1] ) / temp;
  }
  
  b[N-1] = (b[N-1] - lower[N-2]*b[N-2] ) / (diagonal[N-1] - lower[N-2]*upper[N-2] );
  
  // BACKWARD SUBSTITUTION
  x[N-1] = b[N-1];
 
  for (idx i = N-1; i>0 ; i--)
  {
      x[i-1] = (b[i-1] - upper[i-1]*x[i] );
  }

}

void TDMatrix::print() const
{
  /*
  cout << "Tridiagonal matrix size: "<< size << endl;
  for (idx i = 0 ; i < size - 1 ; i++)
    cout << "l,d,u " << i <<" = " << lower[i] << "," << diagonal[i] << "," << upper[i] << endl;
  
  cout << "  d   " << size-1<< " = " << diagonal[size-1] << endl;
  */
  
  for (idx r = 0 ; r < size ; r++)
  {
    for (idx c = 0 ; c < size ; c++)
    {
     if (this->isValidIndex(r,c) )
       printf( "%1.1f\t", this->sparse_get(r,c) );
     else
       printf("%1.1f\t",0.0);
    }
    printf("\n");
  }
    
}


