#include <luincpreconditioner.h>

LUIncPreconditioner::LUIncPreconditioner(const IRCMatrix &A):
    a(A)
{
/*!
 * Creates an incomplere LU preconditioner with zero drop tolerance, LU(0).
 * The implementation is based on the IKJ version of the ILU factorisation 
 * in Saad's book (ALGORITHM 10.3: General ILU Factorization, IKJ Version)
 */
  
  
  
    idx n = A.getNumRows();
    
    // 
    for (idx i = 1 ; i < n ; ++i)
    {
        // FOR ROW k
	for (idx k = 0 ; k < i ; ++k)
        {
            real *aik(NULL);
            if (a.isNonZero(i,k, aik))
            {
                // SCALE L BY DIAGONAL
	        (*aik)/= a.getValue(k,k); // <- TODO : diagonals could be pre-fetched
					  //           in parallel to avoid repeated //           searches.
                
		// MODIFY U. 
		// TODO, OPTIMISE:
		// INDEX j IS OVER A SINGLE ROW, AND THE SEARCH FOR a[i,j] COULD
		// BE AVOIDED BY LOOPING FROM k+1'TH NON-ZERO COLUMN TO TO
		// END OF NON-ZEROS IN ROW i. 
		for (idx j = k+1 ; j < n ; ++j)
                {
                    real *aij(NULL);
                    if ( a.isNonZero(i,j, aij) )
                       (*aij)-= (*aik)*a.getValue(k,j); // COULD NON-ZEROS IN a[k,j] 
                                                        // BE ITERATED OVER TO AVOID SEARCH ?
                }// end for j
            }
        }// end for k
    }// end for i



}

void LUIncPreconditioner::solveMxb(Vector &x, const Vector &b) const
{
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

    // FOR EACH ROW
    for(idx i = 0 ; i < x.getLength() ; ++i)
    {

        // FORM SUM OVER ALL LOWER DIAGONAL VALUES
        real sum(0.0); // ADD IMPLICIT UNIT DIAGONAL ENTRY

        std::vector<IndVal>::const_iterator itr1 = a.nonZeros[i].begin();
        idx col = itr1->ind;
        while ( col < i )
        {
            sum += x[col] * (itr1->val);
            itr1++;
            col = itr1->ind;
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

    idx n = b.getLength(); // NUMBER OF ROWS

    // FOR EACH ROW, COUNT BACKWARDS
    for (idx i = n-1 ; i < n ; --i)
    {
        //real sum(0.0);
        x[i] = b[i];

        std::vector<IndVal>::const_reverse_iterator itr = a.nonZeros[i].rbegin();
        idx col = itr->ind;
        while ( col > i)
        {
            x[i]-=itr->val*x[col];
            itr++;
            col = itr->ind;
        }
        x[i] /= itr->val; // FINALLY DIVISION BY DIAGONAL

    }

}




void LUIncPreconditioner::print() const
{
    a.print();
}
