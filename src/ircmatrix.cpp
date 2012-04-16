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
    return this->cvPairs[i].val;
}