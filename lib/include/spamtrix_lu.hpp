// Copyright (c) 2026 Eero
// SPDX-License-Identifier: MIT
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
/**
 * @brief Non-pivoted sparse LU factorization using a Crout-style algorithm.
 */
    idx numRows;
    FlexiMatrix L;
    FlexiMatrix U;
    LU(){} //
    /**
     * @brief Forward substitution with the lower factor.
     *
     * @param x Solution vector to populate.
     * @param b Right-hand-side vector.
     */
    void forwardSubstitution(Vector& x, const Vector &b) const;
    /**
     * @brief Backward substitution with the upper factor.
     *
     * @param x Solution vector to populate.
     * @param b Right-hand-side vector.
     */
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
    /**
     * @brief Factorize the supplied matrix.
     *
     * @param A Matrix to factorize.
     */
    LU( const IRCMatrix &A);
    /** @brief Destroy the factorization. */
    virtual ~LU();
    /** @brief Print the factorization. */
    void print();
    /**
    * @brief Solve $A x = b$ using forward and backward substitution.
    *
    * @param x Solution vector to populate.
    * @param b Right-hand-side vector.
    */
    void solve(Vector& x, const Vector& b) const;
};
} // end namespace SpaMtrix
#endif // LU_H
