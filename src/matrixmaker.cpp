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
  known, it is added too, otherwise zero is assumed.
*/
    nz.addNonZero(row, col , val);
}


void MatrixMaker::poisson5Point()
{
  /*!
   * Makes a 5 point finite differences (2D) poisson test matrix.
   * The grid spacing is assumed to be unity, resulting in a main
   * diagonal with value 4, and all off-diagonals with values -1.
   *
   * The FD grid is assumed to have an equal number of rows and columns n,
   * where n = sqrt( side length ) of the built matrix A.
   *
   * Example:
   *          To build matrix that corresponds to a 10 x 10 FD grid,
   *          create a MatrixMake object:
   *             1. MatrixMaker mm(100, 100).
   *             2. mm.poisson5Point().
   *             3. IRCMatrix A = mm.getIRCMatrix();
   */

    // MAKE SURE MATRIX IS SQUARE AND NON-ZERO SIZED
    assert(nRows == nCols);
    assert(nRows);

    idx n = sqrt(nRows); // FINITE DIFFERENCES GRID SIDE LENGTH
    // SET MAIN DIAGONALS TO 4
    for (idx i = 0 ; i  < nRows; i++)
        addNonZero(i,i, 4);


    // SET OFF-DIAGONALS
    for (idx i = 0 ; i < nRows ; i++)
    {
        idx row = i / n; // ROW OF i'th NODE
        idx col = i % n; // COLUMN OF i'th NODE

        // RIGHT DERIVATIVE
        if (col < n-1)
            addNonZero(i, i+1, -1);
        // LEFT DERIVATIVE
        if (col>0)
        {
            // cout << "col : " << col << ", i : "<< i <<endl;
            addNonZero(i, i-1, -1);
        }
        // UP DERIVATIVE
        if (row > 0)
        {
            addNonZero(i - n , i, -1);
        }
        if (row < n-1)
        {
            addNonZero(i + n , i , -1);
        }
    }// end off-diagonals



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

