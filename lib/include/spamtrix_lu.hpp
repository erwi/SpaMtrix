#ifndef LU_H
#define LU_H
#include <assert.h>
#include <stdio.h>
#include <omp.h>

#include <spamtrix_setup.hpp>
#include <spamtrix_ircmatrix.hpp>
#include <spamtrix_fleximatrix.hpp>


namespace SpaMtrix{
class LU{
/*!
* A class that performs:\n
*\t1.\t(non-pivoted) LU factorisation of a sparse matrix using the Crout algorithm.\n
*\t2.\tforward-backward bustitution to solve the Ax=b system, where A is the factorised matrix
*/
    idx numRows;
    FlexiMatrix L;
    FlexiMatrix U;
    LU(){} //
    void forwardSubstitution(Vector& x, const Vector &b) const;
    void backwardSubstitution(Vector &x, const Vector &b) const;
    inline void fillFirstColumnL(const IRCMatrix& A, FlexiMatrix &L, const idx &n) {
        for (idx i = 0 ; i < n ; ++i){
            real val;
            if (A.isNonZero(i,0,val)) {
              L.addNonZero(i, 0, val);
            }
        }
    }

    inline void fillFirstRowU(const IRCMatrix &A, FlexiMatrix &U, const idx& n){
        // NORMALISED ROW TO U
        real A00 = A.getValue(0,0);
        for (idx i = 0 ; i < n ; ++i){
            real val;
            if ( A.isNonZero(0,i,val) ){
                U.addNonZero(0,i, val / A00);//L.getValue(0,0) );
            }
        }
    }
public:
    LU( const IRCMatrix &A);
    virtual ~LU();
    void print();
    /**
    * Solves Ax=b using forward/backward substitution. <br>
    * Ax  = b <br>
    * L(Uy) = b <br>
    * Lx = y <br>
    */
    void solve(Vector& x, const Vector& b) const;
};
} // end namespace SpaMtrix
#endif // LU_H
