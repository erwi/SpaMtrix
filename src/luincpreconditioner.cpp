#include <luincpreconditioner.h>
namespace SpaMtrix
{
LUIncPreconditioner::LUIncPreconditioner(const IRCMatrix &A):
    a(A)
{
    /*!
 * Creates an incomplere LU preconditioner with zero drop tolerance, LU(0).
 * The implementation is based on the IKJ version of the ILU factorisation
 * in Saad's book (ALGORITHM 10.3: General ILU Factorization, IKJ Version)
 */

    const idx n = A.getNumRows();

    for (idx i = 1 ; i < n ; ++i)
    {
        std::vector<IndVal>::const_iterator rowEnd = a.nonZeros[i].end();
        std::vector <IndVal> :: iterator k_itr = a.nonZeros[i].begin();
        // k-INDEX LOOP OVER ALL NONZEROS IN ROW i WITH COLUMN INDEX < i
        // (CONSTRUCTS L-PART OF LU)
        while (k_itr->ind < i ) // loop a[i,k]
        {
            const idx k = k_itr->ind;
            k_itr->val /= a.getValue(k,k);
            // j-LOOP OVER NONZEROS IN ROW i WITH COL INDEX > k
            auto j_itr = k_itr + 1;
            while ( j_itr != rowEnd ) // loop a[i,j]
            {
                const idx j = j_itr->ind;
                j_itr->val -= (k_itr->val)*a.getValue(k,j); // COULD THIS SEARCH BE OPTIMISED OUT?
                j_itr ++;
            }// end for/while j
            k_itr ++;
        }// end for/while k
    }// end for i



}

void LUIncPreconditioner::solveMxb(Vector &x, const Vector &b) const
{
#ifdef DEBUG
    assert(x.getLength() == b.getLength() );
#endif
    Vector y( b.getLength() ); // TEMPORARY VECTOR
    forwardSubstitution(y,b);
    backwardSubstitution(x,y);
}

void LUIncPreconditioner::forwardSubstitution(Vector &x, const Vector &b) const
{
    /*!
  Performs forward substitution on the lower diagonal part of the incomplete
  LU factorisation (L).

  L is a lower unit diagonal matrix, stored in the lower part of matrix 'a'.
  The diagonal entries of L are not explicitly stored, since they are all ones.

    x[i] = { b[i] - sum( L[i,j]*x[j] ) } / L[i,j];
    simplifies to
    x[i] = b[i] - sum( L[i,j]*x[j] )

*/
    const idx n = x.getLength();
    std::vector<IndVal>::const_iterator col_itr;
    // FOR EACH ROW
    for(idx i = 0 ; i < n ; ++i)
    {
        // FORM SUM OVER ALL LOWER DIAGONAL VALUES
        real sum(0.0);
        col_itr = a.nonZeros[i].begin();
        while (col_itr->ind < i)
        {
            sum += x[col_itr->ind] * (col_itr->val);
            ++col_itr;
        }
        x[i] = b[i]-sum;
    }
}

void LUIncPreconditioner::backwardSubstitution(Vector &x, const Vector &b) const
{
    /*!
    Performs back-substitution of upper diagonal part of the incomplete
    LU factorisation (U).

    U is the upper diagonal matrix, stored in the upper diagonal part of 'a'.
  */

    const idx n = b.getLength(); // NUMBER OF ROWS
    std::vector<IndVal>::const_reverse_iterator col_itr;
    // FOR EACH ROW, COUNT BACKWARDS
    for (idx i = n-1 ; i < n ; --i)
    {
        x[i] = b[i];
        col_itr = a.nonZeros[i].rbegin();

        while ( col_itr->ind > i)
        {
            x[i]-=col_itr->val*x[col_itr->ind];
            ++col_itr;
        }
        x[i] /= col_itr->val; // FINALLY DIVISION BY DIAGONAL
    }
}




void LUIncPreconditioner::print() const
{
    a.print();
}
} // end namespace