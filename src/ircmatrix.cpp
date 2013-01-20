#include <ircmatrix.h>


namespace SpaMtrix
{
IRCMatrix::IRCMatrix():rows(0), cvPairs(0), nnz(0), numRows(0), numCols(0){/*!Constructs an empty matrix*/}
IRCMatrix::IRCMatrix(  const idx numRows, const idx numCols,
                const idx nnz, 
                idx * const rows, IndVal *const cvPairs):
                rows(rows), cvPairs(cvPairs), nnz(nnz), 
                numRows(numRows) , numCols(numCols){ 
		  /*!Constructs a sparse matrix where sparsity pattern is
		  defined in the arrays 'rows' and 'cvPairs'. The matrix takes ownership of these arrays and releases
		  the memory when destructed (i.e. does not make copies).
		  */
		}
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
        std::cout << "error in " << __func__ << " could not allocate rows" << std::endl;
        exit(1);
    }
    
    cvPairs = new IndVal[nnz];
    if (!cvPairs)
    {
        std::cout << "error in " << __func__ << " could not allocate cvPairs" << std::endl;
        exit(1);
    }
    
    memcpy(rows, m.rows, (nnz+1)*sizeof(idx));
    memcpy(cvPairs, m.cvPairs, nnz*sizeof(IndVal) );
}

IRCMatrix& IRCMatrix::operator=(const IRCMatrix& M)
{
    /*!
      Asigment operator. Sets this matrix equal to matrix M. Reallocates memory if necessary
    */
    if (&M == this)
        return *this;

    numRows = M.numRows;
    numCols = M.numCols;
    // IF OTHER MATRIX SIZE IS DIFFERENT FROM CURRENT, MUST REALLOCATE
    if (nnz != M.nnz)
    {
        // TODO: ADD CHECKS FOR ALLOCATION FAILS
        delete [] rows;
        rows = new idx [nnz+1];

        delete [] cvPairs;
        cvPairs = new IndVal[nnz];

        nnz = M.nnz;
    }
    memcpy(rows, M.rows, (nnz+1)*sizeof(idx));
    memcpy(cvPairs, M.cvPairs, nnz*sizeof(IndVal));

    return *this;
}

IRCMatrix& IRCMatrix::operator=(const real &s)
{
      /*! SETS ALL NONZEROS TO SCALAR s.*/
      for (idx i = 0 ; i < nnz ; ++i)
      {
          cvPairs[i].val = s;
      }
      return *this;
}

const IRCMatrix IRCMatrix::operator*(const real& s) const
{
/*! Matrix scalar multiplication. Returns a scaled version of self.*/  

    idx *rows_n = new idx[nnz+1];
    IndVal *cvPairs_n = new IndVal[nnz];
  
    memcpy(rows_n, rows, (nnz+1)*sizeof(idx) );
    
    for(idx i = 0 ; i < nnz ; ++i)
    {
      cvPairs_n[i] = cvPairs[i];
      cvPairs_n[i].val*=s;
    }

    return IRCMatrix(numRows, numCols, nnz, rows_n, cvPairs_n);
}

idx IRCMatrix::getIndex(const idx row, const idx col) const
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
        if (this->cvPairs[i].ind == col )
            return i;
        ++i;
    }
    // IF REACHED THIS POINT, COLUMN NOT FOUND

    std::cout << "error in " << __func__ << " index " << row <<","<< col <<" not found - bye!" << std::endl;
    exit(1);
}

IRCMatrix::~IRCMatrix()
{
  clear();
}

void IRCMatrix::clear()
{
  /*!
   * Clears all data and deallocates memory.
   */
  delete [] rows;
  delete [] cvPairs;
  rows = NULL;
  cvPairs = NULL;
  numCols = 0;
  numRows = 0;
  nnz = 0;
}

void IRCMatrix::copyFrom(const FlexiMatrix &A)
{
/*!
 * Makes a copy of FlexiMatrix A
 */
    if (nnz > 0 )
      clear();

    nnz = A.calcNumNonZeros();
    numRows = A.numDim1;
    numCols = A.numDim2;
        
    cvPairs = new IndVal[nnz];
    rows = new idx[numRows+1];
    
    idx nzc = 0;    // COUNTER FOR NON-ZEROS
    // COPY DATA
    for (idx r = 0 ; r < numRows ; ++r) // FOR EACH ROW
    {
      rows[r] = nzc; // INDEX TO START OF ROW r
      std::vector<IndVal>::const_iterator citr = A.nonZeros[r].begin();
      for (;citr!=A.nonZeros[r].end(); ++citr)
      {
	cvPairs[nzc] = *citr;
	nzc++;
      }// end column iterator
    }// end for rows
    rows[numRows] = nzc;  
}

void IRCMatrix::sparse_set(const idx row, const idx col, const real val)
{
    idx i = getIndex(row, col);
#ifdef USES_OPENMP
    // #pragma omp critical
#endif
    {
        cvPairs[i].val = val;
    }
    
}

void IRCMatrix::sparse_add(const idx row, const idx col, const real val)
{
    idx i = getIndex(row, col);
#ifdef USES_OPENMP
    //  #pragma omp atomic
#endif
    this->cvPairs[i].val+= val;
}

real IRCMatrix::sparse_get(const idx row, const idx col) const
{
    
    idx i = getIndex(row,col);
    real val = this->cvPairs[i].val;
    return val;
}

real IRCMatrix::getValue(const idx row, const idx col) const
{
    /*!
      Returns value at (row,col). If a non-zero does not exist at
      (row,col), a zero is returned. \n\n
      Note:\n
      Use bool IRCMatrix::isNonZero(row, col, val) instead if it
      is important to know wheteher position (row,col) is a non-zero
      */
    real val;
    isNonZero(row,col,val);
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

    IndVal *begin = &cvPairs[ rows[row]   ]; // pointer to first in row
    IndVal *end   = &cvPairs[ rows[row+1] ]; // pointer to first in row+1


    // GET POINTER TO FIRST ELEMENT IN cvPairs WHOSE INDEX IS
    // LARGER OR EQUAL TO col
    IndVal *itr =
            std::lower_bound(begin,
                             end,
                             col,
                             [](const IndVal &iv1, const idx &colind){return iv1.ind < colind;}
                             );


    if (itr == end) // IF REACHED END OF ROW
        return false;
    else
    if (itr->ind == col)    // IF FOUND
        return true;
    else
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

    val = 0.0;
    IndVal *begin = &cvPairs[ rows[row]  ]; // pointer to first in row
    IndVal *end   = &cvPairs[ rows[row+1]]; // pointer to first in row+1
    // GET POINTER TO FIRST ELEMENT IN cvPairs WHOSE
    // INDEX IS EQUAL OR LARGER THAN col
    IndVal *itr = std::lower_bound(begin,end,col,
                                   [](const IndVal &iv1, const idx& colind){return iv1.ind<colind;}
                                    );
    if (itr == end) // END OF ROW REACHED
    {
        return false;
    }
    else
    if (itr->ind == col) // FOUND IT
    {
        val = itr->val;
        return true;
    }
    else // DID NOT FIND
        return false;

}

void IRCMatrix::operator*=(const real &s)
{
    /*! Modifies matrix by multiplying all values by scalar coefficient s.*/
    for (idx i = 0 ; i < nnz ; i++)
        cvPairs[i].val*=s;
}

Vector IRCMatrix::operator *(const Vector& x) const
{
    /*!
     * MATRIX VECTOR MULTIPLICATION Ax=b.
     * reference to b is returned
     */

#ifdef DEBUG
    if (numCols != x.getLength() )
    {
        std::cout << "error in " << __func__ <<" column count is :" << numCols << 
        ", vector length is :" << x.getLength() << std::endl;
	exit(1);
    }
    //assert(this->getNumCols() == x.getLength() );
#endif

    Vector b( x.getLength() );

    // FOR EACH ROW
//#ifdef USES_OPENMP
//#pragma omp parallel for
    //#endif
    for (idx i = 0 ; i < getNumRows() ; i++)
    {
        // FOR EACH COLUMN
        real r(0);
        const idx row_start = rows[i];
        const idx row_end   = rows[i+1];
        for (idx j = row_start ; j < row_end ; j++)
        {
            const idx col = cvPairs[j].ind;
            r += cvPairs[j].val * x[col];
        }// end for jj
        b[i] = r;

    }//
    return b;
}


// #ifdef DEBUG


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
            real val = this->getValue(row, col );
            printf("%1.3f\t", val );
        }
        printf("\n");
    }
}

} // end namespace SpaMtrix
// #endif
