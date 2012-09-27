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
    
    cvPairs = new IndVal[nnz];
    if (!cvPairs)
    {
        cout << "error in " << __func__ << " could not allocate cvPairs" << endl;
        exit(1);
    }
    
    memcpy(rows, m.rows, (nnz+1)*sizeof(idx));
    memcpy(cvPairs, m.cvPairs, nnz*sizeof(IndVal) );
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
    delete [] rows;
    delete [] cvPairs;
    numCols = 0;
    numRows = 0;
    nnz = 0;

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
      (row,col), a zero is returned.
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
                             [](const IndVal &iv1, const idx &colind){return iv1.val < colind;}
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


Vector IRCMatrix::operator *(const Vector& x) const
{
    /*!
     * MATRIX VECTOR MULTIPLICATION Ax=b.
     * reference to b is returned
     */

#ifdef DEBUG
    assert(this->getNumCols() == x.getLength() );
#endif

    Vector b( x.getLength() );

    // FOR EACH ROW
#ifdef USES_OPENMP
#pragma omp parallel for schedule(static,10000)
#endif
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
            real val = 0.0;
            if (this->isNonZero(row,col) )
                val = this->sparse_get(row, col);

            printf("%1.3f\t", val );
        }
        printf("\n");
    }
    
    
}


// #endif
