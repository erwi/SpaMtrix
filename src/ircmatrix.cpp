#include <ircmatrix.h>

  IRCMatrix::IRCMatrix(const IRCMatrix& m):
  rows(NULL),
  cvPairs(NULL),
  nnz(0),
  numRows(0),numCols(0)
  {
    nnz = m.nnz;
    numRows = m.numRows;
    numCols = m.numCols;
    
    rows = new idx [nnz+1];
    if (!rows)
    {
      cout << "error in " << __func__ << " could not allocate rows" << endl; 
      exit(1);
    }
    
    cvPairs = new ColVal[nnz];
    if (!cvPairs)
    {
      cout << "error in " << __func__ << " could not allocate cvPairs" << endl;
      exit(1);
    }
    
    memcpy(rows, m.rows, (nnz+1)*sizeof(idx));
    memcpy(cvPairs, m.cvPairs, nnz*sizeof(ColVal) );
  }

  IRCMatrix& IRCMatrix::operator=(const IRCMatrix& m)
  {
    // TODO 
      std::cout<< "unimplemented functionality: "<<__func__<< "in file "<< __FILE__<< std::endl;
    return *this;
  }
  
  inline idx IRCMatrix::getIndex(const idx row, const idx col) const
  {
  /*! FINDS INDEX TO ColVal CORRESPONDING INPUT TO ROW AND COLUMN
      TERMINATES WITH ERROR IF NOT FOUND
  */
  #ifdef DEBUG
      assert(row<this->numRows);
      assert(col<this->numCols);
  #endif    
      
      idx i = this->rows[row];
      idx max = this->rows[row+1];
      while ( i < max )
      {
	  if (this->cvPairs[i].col == col )
	      return i;
	  ++i;
      }
      // IF REACHED THIS POINT, COLUMN NOT FOUND
      
      std::cout << "error in " << __func__ << " index " << row <<","<< col <<" not found - bye!" << std::endl;
      exit(1);
  }

  IRCMatrix::~IRCMatrix()
  {
    delete [] rows;
    delete [] cvPairs;
    numCols = 0;
    numRows = 0;
    nnz = 0;
   
  }
  void IRCMatrix::sparse_set(const idx row, const idx col, const real val)
  {
      idx i = getIndex(row, col);
 
  //    #pragma omp atomic
      this->cvPairs[i].val = val;
  }

  void IRCMatrix::sparse_add(const idx row, const idx col, const real val)
  {
      idx i = getIndex(row, col);
  //    #pragma omp atomic
      this->cvPairs[i].val+= val;
  }

  real IRCMatrix::sparse_get(const idx row, const idx col) const
  {
    
      idx i = getIndex(row,col);
      real val = this->cvPairs[i].val;
      return val;
  }

bool IRCMatrix::isNonZero(const idx row, const idx col) const
{
/*! Returns true if this matrix contains storage at location row, col, otherwise false
*/
#ifdef DEBUG
  assert(row < this->numRows);
  assert(col < this->numCols);
#endif
  idx i = this->rows[row];
  idx max = this->rows[row+1];
  while ( i < max )
  {
    if (this->cvPairs[i].col == col )
      return true;
    ++i;
  }
  // IF REACHED THIS POINT, COLUMN NOT FOUND
  return false;
  }

bool IRCMatrix::isNonZero(const idx row, const idx col, real &val) const
{
/*!
 Returns true if this matrix contains storage at location row, col, otherwise false.
 The actual value of location row,col is returned in val
*/
#ifdef DEBUG
    assert(row < this->numRows);
    assert(col < this->numCols);
#endif
    idx i = this->rows[row];
    idx max = this->rows[row+1];

    while (i<max)
    {
        if (this->cvPairs[i].col == col)
        {
            val = cvPairs[i].val;
            return true;
        }
        ++i;
    }

    val = 0.0;
    return false;
}


  #ifdef DEBUG


  void IRCMatrix::spy() const
  {
    std::cout << std::endl;
    std::cout << "IRCMatrix size = " << this->numRows << " , " << this->numCols << " nnz = "<< this->nnz << std::endl;
    for (idx row = 0 ; row < this->numRows ; row++)
    {
      for (idx col = 0 ; col < this->numCols ; col++)
      {
	char marker;
	if ( this->isNonZero(row,col) )
	  marker = '#';
	else
	  marker = '.';
	printf("%c", marker );
	
      }
      printf("\n");
      
    }
  }
  
  void IRCMatrix::print() const
  {
    std::cout << std::endl;
    std::cout << "IRCMatrix size = " << this->numRows << " , " << this->numCols << " nnz = "<< this->nnz << std::endl;
    for (idx row = 0 ; row < this->numRows ; row++)
    {
      for (idx col = 0 ; col < this->numCols ; col++)
      {
	  real val = 0.0;
	if (this->isNonZero(row,col) )
	  val = this->sparse_get(row, col);
 	
	printf("%1.3f\t", val ); 
      }
      printf("\n");
    }
    
    
  }
  
  
  #endif
