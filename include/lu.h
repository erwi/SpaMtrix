#ifndef LU_H
#define LU_H
#include <assert.h>
#include <setup.h>
#include <ircmatrix.h>
#include <fleximatrix.h>
#include <stdio.h>
#include <omp.h>
/*!
  A class that performs (non-pivoted) LU factorisation of a sparse matrix using the Crout algorithm.
  */
namespace SpaMtrix
{
class LU
{
    idx numRows;

    FlexiMatrix L;
    FlexiMatrix U;

    LU(){} //

    void forwardSubstitution(Vector& x, const Vector &b) const;
    void backwardSubstitution(Vector &x, const Vector &b) const;

    inline void fillFirstColumnL(const IRCMatrix& A, FlexiMatrix &L, const idx &n)
    {
        for (idx i = 0 ; i < n ; ++i)
        {
            real val;
            if ( A.isNonZero(i,0,val) )
                L.addNonZero(i,0,val);
        }
    }

    inline void fillFirstRowU(const IRCMatrix &A, FlexiMatrix &U, const idx& n)
    {
        // NORMALISED ROW TO U
        real A00 = A.getValue(0,0);
        for (idx i = 0 ; i < n ; ++i)
        {
            real val;
            if ( A.isNonZero(0,i,val) )
            {
                U.addNonZero(0,i, val / A00);//L.getValue(0,0) );
            }
        }
    }


public:
    LU( const IRCMatrix &A);
    void print();
    void solve(Vector& x, const Vector& b) const; // SOLVES Ax = b USING FORWARD/BACKWARD SUBSTITUTION
};
} // end namespace SpaMtrix
#endif // LU_H
