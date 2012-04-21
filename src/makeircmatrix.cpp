#include <ircmatrix.h>
#include <algorithm>
#include <iostream>
#include <vector>
#include <cstring>
using std::cout;
using std::endl;
IRCMatrix& makeIRCMatrix(const idx* t,
                         const idx npt,
                         const idx nt)
{
/*! 
    Makes an IRCMatrix based on the connectivity defined in the array t
*/

// DETEMINE REQUIRED MATRIX SIZE. NUMBER OF ROWS/COLUMNS IS 
// LARGEST VALUE FOUND IN t
    idx numRows = *std::max_element(t , t+(npt * nt ) ) + 1;
    
    // ALLOCATE VECTOR OF VECTORS OF COLUMN INDEXES
    std::vector < std::vector < idx > > col_idxs(numRows);
    for (idx i = 0 ; i < numRows ; ++i) // INITIALISE TO EMPTY
    {
        std::vector <idx> empty;
        col_idxs.push_back(empty);
    }

    // MAKE COLUMN INDEX ARRAYS
    for (idx i = 0 ; i < nt ; i++)// FOR EACH ELEMENT i
    {
        idx offset = i*npt; // INDEX TO START OF ELEMENT
        for (idx j = 0 ; j < npt ; j++ )// FOR EACH NODE j IN ELEMENT i
        {
            idx n1 = t[offset + j];
            for (idx k = 0 ; k < npt ; k++ )// FOR EACH NODE k IN ELEMENT i
            {
                idx n2 = t[offset + k];
                col_idxs[n1].push_back(n2); // ADD COLUMN INDEXES
                col_idxs[n2].push_back(n1);
            }
        }        
    }// END FOR i
    
    // SORT AND UNIQUEFY COLUMN INDEX ARRAYS
    #pragma omp parallel for
    for (idx i = 0 ; i < numRows ; i++ )
    {
        std::sort( col_idxs[i].begin() , col_idxs[i].end() );
        std::vector<idx> :: iterator itr;
        
        itr = std::unique( col_idxs[i].begin(), col_idxs[i].end() );
        col_idxs[i].resize(itr - col_idxs[i].begin() );
    }
    
    // FIND TOTAL NUMBER OF NON-ZEROS
    idx nnz = 0;
    for (idx i = 0 ; i < numRows ; i++ )
    {
        nnz += (idx) col_idxs[i].size();
    }
        
    // ALLOCATE MATRIX INDEX ARRAYS
    idx *rows = new idx[numRows + 1];	// ROW COUNTER 
    ColVal* cvPairs = new ColVal[nnz];	// COLUMN INDEXES (AND DATA)
    memset(cvPairs, 0 , nnz*sizeof( ColVal ) );
    
    // FILL IN COLUMN INDEX DATA AND ROW COUNTER ARRAY 
    idx c = 0;
    for (idx i = 0 ; i < numRows ; i++ ) // FOR EACH ROW
    {
        rows[i] = c;	// SET INDEX TO START OF ROW i
        idx numC = col_idxs[i].size(); // NNZs IN ROW i
        for (idx j = 0 ; j < numC ; j++ )
        {
           // COPY VALUE FROM jth COLUMN INDEX TO cvPair 
	  (cvPairs[c]).col = *(col_idxs[i].begin()+j); 
           c++;
        }
    }
    rows[numRows] = c; // LAST NNZ+1
    
        
    // FINALLY CREATE MATRIX
    IRCMatrix* ircm = new IRCMatrix(numRows, numRows , nnz, rows, cvPairs);
    return *ircm;

}