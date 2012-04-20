#include <ircmatrix.h>

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


  void IRCMatrix::sparse_set(const idx row, const idx col, const real val)
  {
      idx i = getIndex(row, col);
      #pragma omp atomic
      this->cvPairs[i].val = val;
  }

  void IRCMatrix::sparse_add(const idx row, const idx col, const real val)
  {
      idx i = getIndex(row, col);
      #pragma omp atomic
      this->cvPairs[i].val+= val;
  }

  real IRCMatrix::sparse_get(const idx row, const idx col) const
  {
    
      idx i = getIndex(row,col);
      real val = this->cvPairs[i].val;
      return val;
  }


  #ifdef DEBUG

  bool IRCMatrix::isNonZero(const idx row, const idx col) const
  {
  /*! Returns true if this matrix contains storage at location row, col, otherwise false
  */
    assert(row < this->numRows );
    assert(col < this->numCols );
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
