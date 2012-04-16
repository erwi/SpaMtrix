#ifndef IRCMATRIX
#define IRCMATRIX
#include <iostream>
#include <assert.h>
typedef unsigned int idx;
typedef double real;
struct ColVal{
	idx col;	// COLUMN INDEX
	real val;	// VALUE
};



// Interleaved Row Compressed Matrix
class IRCMatrix
{
    idx* rows;			// ROW COUNTER
    ColVal* cvPairs;	// COLUMN-VALUE PAIRS
    idx nnz;			// NUMBER OF NON-ZEROS
    idx numRows, numCols;
    
    idx getIndex(const idx row, const idx col) const;
public:
    
    IRCMatrix():
        rows(0), cvPairs(0), nnz(0), numRows(0), numCols(0) {}
    IRCMatrix(  const idx numRows, const idx numCols,
                const idx nnz, 
                idx * const rows, ColVal *const cvPairs):
                rows(rows), cvPairs(cvPairs), nnz(nnz), 
                numRows(numRows) , numCols(numCols){}
    
    //================================================
    idx getnnz()const {return nnz;}
    idx getNumRows()const {return numRows;}
    idx getNumCols()const {return numCols;}
    void sparse_set(const idx row, const idx col , const real val );
    void sparse_add(const idx row, const idx col , const real val );
    real sparse_get(const idx row, const idx col ) const;
    
    
};

// MAKES SIMPLE IRCMatrix OUT OF AN ARRAY OF CONNECTED INDEXES
IRCMatrix& makeIRCMatrix(const idx* t, const idx npt, const idx nt);

#endif
