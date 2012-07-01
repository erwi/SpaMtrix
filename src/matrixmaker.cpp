#include <matrixmaker.h>

MatrixMaker::MatrixMaker():
nRows(0),
nCols(0)
{
  
}

void MatrixMaker::setMatrixSize(const idx nRows, const idx nCols)
{
  /*!
   * Allocates space for the lists of non-zeros. Use this before starting
   * entering non-zeros
   */
  
  // CLEAR OLD DATA, IF ANY
  if (!nonZeros.empty() )
    nonZeros.clear();
  
  this->nRows = nRows;
  this->nCols = nCols;
  
  // ALLOCATE NONZEROS AND AUTOMATICALLY ADD DIAGONAL COMPONENT
  for (idx i = 0 ; i < nRows ; i++)
  {
    std::list<idx> l;
    l.push_back(i);
    nonZeros.push_back(l);
  }
  
  
}

void MatrixMaker::addNonZero(const idx row, const idx col)
{
  /*!
   * Adds storage for a non-zero at location row,col.
   */
  
  if ( ( row >= this->nRows) || (col >= this->nCols) )
  {
    cout << "error in " << __func__ << "cannot add nonzero at:" << row << "," << col << endl;
    cout << "Matrix size was earlier defined as : "<< nRows << "," << nCols << endl;
    exit(1);
  }
  
  this->nonZeros[row].push_back(col);
   
}

void MatrixMaker::clear()
{
  nRows = 0;
  nCols = 0;
  nonZeros.clear();
}

IRCMatrix MatrixMaker::getIRCMatrix()
{
  /*!
   * Converts sparsity pattern to 0-initialised IRCMatrix
   */
  // SORT AND UNIQUEFY COLUMN INDEXES
  for (idx i = 0 ; i < nRows ; i++)
  {
    // SORT
    nonZeros[i].sort();
    // REORDER
    std::list<idx> :: iterator itr;
    itr = std::unique(nonZeros[i].begin() , nonZeros[i].end() );
    // REMOVE REPEATED
    nonZeros[i].erase(itr, nonZeros[i].end() );
  }
  
  // FIND TOTAL NUMBER OF NON-ZEROS
  idx nnz = 0;
  for (idx i = 0 ; i < nRows ; i++)
  {
    nnz += (idx) nonZeros[i].size();
  }
  
  // ALLOCATE MEMORY FOR FINAL MATRIX ARRAYS
  idx* rows = new idx[nRows + 1]; 	// ROW COUNTER
  ColVal* cvPairs = new ColVal[nnz]; 	// COLUMN INDEXES AND VALUES
  memset(cvPairs, 0 , nnz*sizeof( ColVal ) );
  
  // FILL IN COLUMN AND ROW INDEX ARRAYS
  idx cnt = 0; // counter
  for (idx r = 0 ; r < nRows ; r++)
  {
    rows[r] = cnt; // INDEX TO START OF ROW r
 
    // FOR EACH NON-ZERO IN ROW r
    list <idx> :: iterator citr;
    for (citr = nonZeros[r].begin() ; citr != nonZeros[r].end(); citr++) 
    {
      // COPY INDEX FROM jTH COLUMN INDEX TO cvPairs
      cvPairs[cnt].col = *citr;
      cnt++;
    }// end for cols
  }// end for rows
  rows[nRows] = cnt; // LAST = NNZ+1
  
  // CREATE MATRIX - FINGERS CROSSED FOR RETURN VALUE OPTIMIZATION
  // USE C++11 MOVE SEMANTICS HERE?
  return IRCMatrix(nRows, nCols, nnz, rows, cvPairs);
 }


