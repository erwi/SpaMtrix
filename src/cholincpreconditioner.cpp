#include <cholincpreconditioner.h>

CholIncPreconditioner::CholIncPreconditioner(const IRCMatrix& A):
    LD(0)
{
    /*!
 * Incoplete cholesky matrix constructor.
 */
    A.spy();
    makeLowerDiagonal(A);
    LD->spy();


}

void CholIncPreconditioner:: makeLowerDiagonal(const IRCMatrix &A)
{
    /*!
  Creates empty sparse matrix with same sparsity pattern as lower half of A
*/
    // FIRST CREATE SPERSITY PATTERN FOR LOWER TRIANGULAT MATRIX OF A
    idx numRows = A.getNumCols();
    idx numCols = A.getNumCols();
    idx nnzA    = A.getnnz();


    // NON-ZEROS FOR THIS WILL  INCLUDE LOWER HALF AND DIAGONAL
    // ENTRIES OF MATRIX A
    idx nnz_diag = numRows;
    idx nnz_lower = ( nnzA - nnz_diag ); // this should be divisible by 2 (always?)
#ifdef DEBUG
    if (nnz_lower % 2 )
    {
        std::cout << "error in " << __func__ << "matrix sparsity doesn't seem to be symmetric - bye!" << std::endl;
        exit(1);
    }
#endif
    idx nnz = ( nnz_lower / 2 ) + nnz_diag; // FINAL NUMBER OF NONZEROS

    // ALLOCATE MEMORY FOR ROW-COUNTER AND COLUMN INDEX/VALUE PAIRS
    idx* rows = new idx[numRows + 1];
    IndVal* cvPairs = new IndVal[nnz];
    memset(cvPairs, 0, nnz*sizeof(IndVal) );

    // FILL IN COLUMN AND ROW INDEX ARRAYS
    // LOOP OVER EACH ROW OF A, SELECTING ONLY COLUMNS THAT ARE BELOW DIAGONAL
    idx cnt = 0; // COUNTER OF NONZEROS
    for (idx row = 0 ; row < numRows ; row++)
    {
        rows[row] = cnt;

        idx rstart = A.rows[row];   // INDEX TO STARTS OF THIS AND NEXT ROWS IN A
        idx rend = A.rows[row+1];

        for (idx i = rstart ; i < rend ; i++) // FOR EACH NONZERO COLUMN IN THIS ROW
        {
            IndVal *cv = &A.cvPairs[i];
            if ( cv->ind > row ) // IF UPPER DIAGONAL COLUMN, THIS ROW IS DONE, START NEXT
            {
                break;
            }
            else //  LOWER DIAGONAL -> COPY COLUMN INDEX AND INCREMENT COUNTER
            {
                cvPairs[cnt].ind = cv->ind;
                cnt++;
            }

        }
    }
    rows[numRows] = cnt;

    LD = new IRCMatrix(numRows, numCols, nnz, rows, cvPairs);

    if (!LD)
    {
        std::cout << "error in " << __func__ << "received null pointer - bye!" << std::endl;
        exit(1);
    }
}


void CholIncPreconditioner::solveMxb(Vector &x, const Vector &b) const
{

}
