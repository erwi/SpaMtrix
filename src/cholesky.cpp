#include <cholesky.h>

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
        L.nonZeros.push_back(vector<IndVal>() ); // ADD NEW EMPTY ROW
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

                s = sqrt(A.sparse_get(r,r) - s );

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
                real s = 0;
                for (idx k = 0 ; k < c; k++)
                  s+= L.getValue(r,k)*L.getValue(c,k);

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
        std::vector<IndVal>::const_iterator itr1 = L.nonZeros[i].begin();
        idx col = itr1->ind;
        // WHILE LOWER DIAGONAL NONZEROS ONLY
        // UPPER DIAGONAL VALUES TAKEN CARE OF IN back-substitution
        while (col < i)
        {
            sum+= x[col] * ( itr1->val );
            itr1++;
            col = itr1->ind;
        }
        x[i] = ( b[i]-sum ) / itr1->val;
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
        real sum(0.0);
        x[i] = b[i];
        // PERFORM SUM OF L'x, WHERE L' IS TRANSPOSE OF L
        // THIS USES MATRIX SEARCH AND SHOULD BE OPTIMISED
        for ( idx j = n-1 ; j > i ; --j)
        {
            x[i]-= L.getValue(j,i)*x[j];
        }
        x[i] /= L.getValue(i,i);
    }// end for each row
}
