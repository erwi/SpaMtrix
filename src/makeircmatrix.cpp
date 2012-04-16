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
cout << "npt = "<< npt << " nt = " << nt << endl;
// DETEMINE REQUIRED MATRIX SIZE. NUMBER OF ROWS/COLUMNS IS 
// LARGEST VALUE FOUND IN t
    idx numRows = *std::max_element(t , t+(npt * nt ) ) + 1;
    std::cout<< "numRows " << numRows << std::endl;

    std::vector < std::vector < idx > > col_idxs(numRows);
    for (idx i = 0 ; i < numRows ; ++i) // INITIALISE
    {
        std::vector <idx> empty;
        col_idxs.push_back(empty);
    }

    // MAKE COLUMN INDEX ARRAYS
    for (idx i = 0 ; i < nt ; i++)
    {
        idx offset = i*npt;
        cout << "i = " << i << " offset = "<< offset << endl;
        for (idx j = 0 ; j < npt ; j++ )
        {
            idx n1 = t[offset + j];
            for (idx k = 0 ; k < npt ; k++ )
            {
                idx n2 = t[offset + k];
                
                cout << "j,k = "<< j << "," << k <<" n1,n2 = " << n1 <<","<<n2<<endl;
                col_idxs[n1].push_back(n2);
                col_idxs[n2].push_back(n1);
            }
        }        
    }
    
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
    cout << "nnz = "<< nnz << endl;
    
    // ALLOCATE MATRIX INDEX ARRAYS
    idx *rows = new idx[numRows + 1];
    ColVal* cvPairs = new ColVal[nnz];
    memset(cvPairs, 0 , nnz*sizeof( ColVal ) );
    // FILL IN CONNECTIVITY DATA
    idx c = 0;
    for (idx i = 0 ; i < numRows ; i++ ) // FOR EACH ROW
    {
        rows[i] = c;
        idx numC = col_idxs[i].size();
        for (idx j = 0 ; j < numC ; j++ )
        {
            (cvPairs[c]).col = *(col_idxs[i].begin()+j);
            c++;
        }
    }
    rows[numRows] = c;
    
    // DEBUG PRINTOUT
    for (idx i = 0; i < numRows+1 ; i++)
    {
        cout << "rows "<< i << " = "<< rows[i] << endl;
    }
    
    for (idx i = 0 ; i < nnz ; i++ )
    {
        cout << cvPairs[i].col <<" ";
    }
    
    
    
    // 
    IRCMatrix* ircm = new IRCMatrix(numRows, numRows , nnz, rows, cvPairs);
    return *ircm;

}