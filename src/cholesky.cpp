#include <cholesky.h>

Cholesky::Cholesky(const IRCMatrix &A):
L(NULL)
{
/*!
  Creates a cholesky decomposition lower diagonal of matrix A
  */
    idx numRows(A.getNumRows());

    std::vector < std::vector < ColVal > > nonZeros( numRows );


    for (idx r = 0 ; r < numRows ; r++)
        for (idx c = 0 ; c <= r ; c++)
        {
            if (r==c) // DIAGONAL
            {
                real s = 0;
                // SUM OVER ALL VALUES SQUARED TO LEFT ON THIS ROW
                std::vector< ColVal>::iterator itr1 = nonZeros[r].begin();
                std::vector< ColVal>::iterator itr2 = nonZeros[r].end();
                for( ; itr1 != itr2 ; itr1++)
                    s+= itr1->val * itr1->val;

                s = sqrt(A.sparse_get(r,r) - s );

                #ifdef DEBUG
                  if (s <= 0.0)
                  {
                    std::cout << "error in " <<__func__ << "matrix A is not positive definite - bye!" << std::endl;
                    exit(1);
                  }
                #endif
                  ColVal cv = {c,s}; // std::move could be used here?
                  nonZeros[r].push_back( cv );


            }
            else // OFF-DIAGONAL
            {
                real s = 0;
                for (idx k = 0 ; k <= c; k++)
                  s+= getValue(r,k, nonZeros)*getValue(c,k, nonZeros);

                real a(0.0);
                A.isNonZero(r,c,a); // GETS VALUE IN A AT r,c. INCLUDING ZERO VALUES


                s =  ( a - s );

                if ( s != 0.0 ) // IF NOT ZERO, ADD TERM TO MATRIX
                {
                    real last= nonZeros[c].back().val;
                          //[c].val;
                  s /= last;
                  ColVal cv = {c,s};
                  nonZeros[r].push_back(cv);
                }
            }
        }




// CONVERT ROWS OF COLUMN-VALUE PAIRS TO AN IRCMatrix
  idx nnz(0);
  for (idx i = 0 ; i < numRows ; i++)
    nnz+= nonZeros[i].size();

  std::cout<< "nnz : " << nnz << std::endl;

  // ALLOCATE MEMORY FOR FINAL MATRIX ARRAYS
  idx* rows = new idx[numRows + 1]; 	// ROW COUNTER
  ColVal* cvPairs = new ColVal[nnz]; 	// COLUMN INDEXES AND VALUES
  assert(rows);
  assert(cvPairs);

  memset(cvPairs, 0 , nnz*sizeof( ColVal ) );

  // FILL IN COLUMN AND ROW INDEX ARRAYS
  idx cnt = 0; // counter
  for (idx r = 0 ; r < numRows ; r++)
  {
    rows[r] = cnt; // INDEX TO START OF ROW r

    // FOR EACH NON-ZERO IN ROW r
    std::vector<ColVal> :: iterator itr1 = nonZeros[r].begin();
    std::vector<ColVal> :: iterator itr2 = nonZeros[r].end();
    for (; itr1!=itr2 ; itr1++)
    {
      // COPY INDEX AND VALUE
      cvPairs[cnt] = *itr1;
      cnt++;
    }
  }// end for rows
  rows[numRows] = cnt; // LAST = NNZ+1

  // CREATE MATRIX - FINGERS CROSSED FOR RETURN VALUE OPTIMIZATION
  // USE C++11 MOVE SEMANTICS HERE?
  L = new IRCMatrix(numRows, numRows, nnz, rows, cvPairs);
  L->print();

}


real Cholesky::getValue(const idx r,
                        const idx c,
                        const std::vector< std::vector <ColVal> >&  L) const
{
/*!
  Returns value of row r, column c of matrix represented by rows of non-zero ColVal
  pairs, L
  */

    // SEARCH THROUGH ROW r FOR AN NON-ZERO ENTRY IN COLUMN c
    // COLUMNS ARE SORTED, SO A BINOMIAL SEARCH COULD BE FASTER HERE
    std::vector<ColVal>::const_iterator itr1 = L[r].begin();
    std::vector<ColVal>::const_iterator itr2 = L[r].end();

    for ( ; itr1 != itr2 ; itr1++)
      if (itr1->col == c)   // IF SOUGHT COLUMN IS FOUND
        return itr1->val;

    // IF NO ENTRY FOUND, RETURN ZERO
    return 0.0;
}
