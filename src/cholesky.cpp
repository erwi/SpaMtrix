#include <cholesky.h>

namespace SpaMtrix
{
Cholesky::Cholesky(const IRCMatrix &A)
{
/*!
  Creates a cholesky decomposition lower diagonal of matrix A
  */
    L.numDim1 = A.getNumRows();
    L.numDim2 = A.getNumCols();

    // FOR EACH ROW
    for (idx r = 0 ; r < A.getNumRows() ; r++)
    {


        L.nonZeros.push_back(std::vector<IndVal>() ); // ADD NEW EMPTY ROW
        real Arr = A.sparse_get(r,r);
        // FOR EACH COLUMN, LOWERD DIAGONAL ONLY, i.e. c < r
        for (idx c = 0 ; c <= r ; c++)
        {
            if (r==c) // DIAGONAL
            {
                real s = 0;
                // SUM OVER ALL VALUES SQUARED TO LEFT ON THIS ROW
                std::vector< IndVal>::iterator itr1 = L.nonZeros[r].begin();
                std::vector< IndVal>::iterator itr2 = L.nonZeros[r].end();
                for( ; itr1 != itr2 ; itr1++)
                    s+= itr1->val * itr1->val;

                s = sqrt(Arr - s );

                #ifdef DEBUG
                  if (s <= 0.0)
                  {
                    std::cout << "error in " <<__func__ << "matrix A is not positive definite - bye!" << std::endl;
                    exit(1);
                  }
                #endif

                  L.nonZeros[r].push_back(IndVal(c,s) );// std::move could be used here?
            }
            else // OFF-DIAGONAL
            {
                real s(0.0);
#ifdef USES_OPENMP
#pragma omp parallel for reduction (+:s)
#endif
                for (idx k = 0 ; k < c; k++)
                {
                    real Lrk = L.getValue(r,k);
                    if ( Lrk != 0.0)
                        s+= Lrk*L.getValue(c,k);
                }


                /*
                const vector<IndVal> &cvPairs = L.nonZeros[r]; // ref to non-zeros in row r
                idx n = (idx) cvPairs.size() - 1;
                for (idx k = 0 ; k < n ; ++k)
                {
                    idx     col = cvPairs[k].ind;
                    real    val = cvPairs[k].val;
                    s+= val*L.getValue(c,col);
                }
                */

                real a(0.0);
                A.isNonZero(r,c,a); // GETS VALUE IN A AT r,c. INCLUDING ZERO VALUES
                s =  ( a - s );

                if ( s != 0.0 ) // IF NOT ZERO, ADD TERM TO MATRIX
                {
                    real last= L.nonZeros[c].back().val;
                    s /= last;
                    L.nonZeros[r].push_back( IndVal(c,s) );
                }
            }
        }//end for columns
    }// end for rows

   // this->makeUpper();
}


void Cholesky::print() const
{
    /*!
      Prints the colesky decmoposed lower triangular matrix to stdout
      */
    L.print();
}

void Cholesky::solve(Vector &x, const Vector &b) const
{
/*!
  Solves Ax=b using forward/backward substitution
    Ax = b
    L'(Ly) = b
    L'x = y
*/
    Vector y( b.getLength() ); // TEMPORARY VECTOR
    forwardSubstitution(y,b);
    backwardSubstitution(x,y);
}

void Cholesky::forwardSubstitution(Vector &x, const Vector &b) const
{
/*!
     x[i] =  ( b[i] - sum( A[i,j]x[j] ) ) / A[i,i]
 */
    // FOR EACH ROW, STARTING FROM FIRST
    for (idx i = 0 ; i < x.getLength() ; ++i)
    {

        // FORM SUM OVER ALL LOWER DIAGONAL MATRIX VALUES
        real sum(0.0);
        // REFERENCE TO ALL NON_ZEROS IN ROW i
        const std::vector<IndVal> &cvPairs = L.nonZeros[i];

        idx n = (idx) cvPairs.size()-1;
        // SUM OVER ALL BELOW DIAGONAL. i.e. DON NOT INCLUDE DIAGONAL
#ifdef USES_OPENMP	
#pragma omp parallel for reduction(+:sum)
#endif
        for (idx j = 0 ; j < n ; ++j)
        {
            const idx col = cvPairs[j].ind;
            const real val = cvPairs[j].val;
            sum += x[col] * val;
        }

        // DIVIDE BY DIAGONAL VALUE
        x[i] = (b[i]-sum) / cvPairs[n].val;

    }
}

void Cholesky::backwardSubstitution(Vector &x, const Vector &b) const
{
/*!
  Performs back substitution of transposed lower triangular matrix.

  x[i] = ( b[i] - sum(L'x[i+1:end]) ) / L(i,i)
*/

    idx n = b.getLength();
    // FOR EACH ROW, STARTING FROM LAST, COUNTING BACKWARDS
    // NOTE: i < n FOR UNSIGNED INTS IS EQUIVALENT TO
    // i >= 0 FOR SIGNED INTS
    for (idx i = n-1 ; i < n ; --i)
    {
        //real sum(0.0);
        //x[i] = b[i];
        // PERFORM SUM OF L'x, WHERE L' IS TRANSPOSE OF L
        // THIS USES MATRIX SEARCH AND SHOULD BE OPTIMISED
        real xi = b[i];
#ifdef USES_OPENMP	
#pragma omp parallel for reduction(-:xi)
#endif
        //for ( idx j = n-1 ; j > i ; --j)
        for (idx j = i+1 ; j < n ; ++j )
        {
            xi-=L.getValue(j,i)*x[j];
        }
        x[i] = xi/L.getValue(i,i);
    }// end for each row
}
} // end namespace SpaMtrix