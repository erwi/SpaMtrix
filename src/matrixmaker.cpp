#include <matrixmaker.h>

/*
MatrixMaker::MatrixMaker():
nRows(0),
nCols(0)
{
  
}
*/
MatrixMaker::MatrixMaker(const idx nRows, const idx nCols):
nRows(nRows),
nCols(nCols),
nz(nRows,nCols)
{

}



void MatrixMaker::addNonZero(const idx row, const idx col, const real val)
{
/*!
  Adds storage for a non-zero at location row, col. If value of non-zero is
  known, it can be added too.
*/

    nz.addNonZero(row, col , val);

}


IRCMatrix MatrixMaker::getIRCMatrix()
{
  /*!
   * Converts sparsity pattern to 0-initialised IRCMatrix
   */

  // FIND TOTAL NUMBER OF NON-ZEROS
  idx nnz = 0;
  for (idx i = 0 ; i < nRows ; i++)
  {
    nnz += (idx) nz.nonZeros[i].size();
  }

  // ALLOCATE MEMORY FOR FINAL MATRIX ARRAYS
  idx* rows = new idx[nRows + 1]; 	    // ROW COUNTER
  IndVal* cvPairs = new IndVal[nnz];    // COLUMN INDEXES AND VALUES
  memset(cvPairs, 0 , nnz*sizeof( IndVal ) );
  //*
  // FILL IN COLUMN AND ROW INDEX ARRAYS
  idx cnt = 0; // counter
  for (idx r = 0 ; r < nRows ; r++)
  {
    rows[r] = cnt; // INDEX TO START OF ROW r
 
    // FOR EACH NON-ZERO IN ROW r
    std::vector <IndVal> :: iterator itr;
    for (itr = nz.nonZeros[r].begin() ; itr != nz.nonZeros[r].end(); itr++)
    {
      // COPY INDEX FROM jTH COLUMN INDEX TO cvPairs
      cvPairs[cnt] = *itr;
      cnt++;
    }// end for cols
  }// end for rows
  rows[nRows] = cnt; // LAST = NNZ+1
  
  // CREATE MATRIX - FINGERS CROSSED FOR RETURN VALUE OPTIMIZATION
  // USE C++11 MOVE SEMANTICS HERE?
  return IRCMatrix(nRows, nCols, nnz, rows, cvPairs);
 }


