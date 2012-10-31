#ifndef SPAMTRIX_BLAS_H
#define SPAMTRIX_BLAS_H

#include <setup.h>
#include <ircmatrix.h>
#include <tdmatrix.h>
#include <vector.h>
#include <assert.h>
#include <iostream>
#include <stdlib.h>
#include <omp.h>


namespace SpaMtrix
{

// MATRIX VECTOR MULTIPLICATION b = Ax
inline void multiply(const IRCMatrix& A,
                     const Vector& x,
                     Vector& b)
{
   /*!
   * Matrix-Vector multiplication Ax=b.
   * \n
   * In python wrappers this function is renamed to IRCMatMul.
   */

#ifdef DEBUG
    assert(A.getNumCols() == x.getLength() );
    assert(b.getLength() == x.getLength() );
#endif
    // FOR EACH ROW
    const idx n = A.getNumRows();
#ifdef USES_OPENMP
#pragma omp parallel for
#endif
    for (idx i = 0 ; i < n ; ++i)
    {
        // FOR EACH COLUMN
        register real r(0);
        const idx row_start = A.rows[i];
        const idx row_end   = A.rows[i+1];
        for (idx j = row_start ; j < row_end ; ++j)
        {
            const idx col = A.cvPairs[j].ind;
            r += A.cvPairs[j].val * x[col];
        }//j
        b[i] = r;

    }//
}

inline real multiply_dot(const IRCMatrix& A,
                     const Vector& x,
                     Vector& b)
{
    /*!
   * MATRIX VECTOR MULTIPLICATION Ax=b
     also returns dot product of x an b. Used in IterativeSolvers
   */

#ifdef DEBUG
    assert(A.getNumCols() == x.getLength() );
    assert(b.getLength() == x.getLength() );
#endif
    // FOR EACH ROW
    const idx n = A.getNumRows();

    real dp(0);
#ifdef USES_OPENMP
#pragma omp parallel for reduction(+:dp)
#endif
    for (idx i = 0 ; i < n ; ++i)
    {
        // FOR EACH COLUMN
        register real r(0);
        const idx row_start = A.rows[i];
        const idx row_end   = A.rows[i+1];
        for (idx j = row_start ; j < row_end ; ++j)
        {
            const idx col = A.cvPairs[j].ind;
            r += A.cvPairs[j].val * x[col];
        }//j
        b[i] = r;
        dp+= b[i]*x[i];
    }//
    return dp;
}


inline void multiply(const TDMatrix& A,
                     const Vector& x,
                     Vector& b)
{
    /*!
   * Matrix vector muliplication Ax = b
   * where A is a tridiagonal matrix \n \n
   * In python this function is renamed to TDMatMul. 
   */

#ifdef DEBUG
    assert(A.size == x.getLength() );
#endif
    idx n = b.getLength();
    for (idx i = 0 ; i < n ; i++)
    {
        b[i] = A.diagonal[i]*x[i];
        if (i > 0 )
           b[i]   += A.lower[i-1] * x[i-1];
        
        if (i < (n-1) )
           b[i] += A.upper[i] * x[i+1];
    }
}
// ======================================
// BLAS LEVEL 1 FUNCTIONS
// ======================================

inline void scale(const real a, Vector& v)
{
    /*!
   *  v*= a
   */
    for (idx i = 0 ; i < v.getLength() ; ++i)
        v[i]*=a;
}


// VECTOR DOT PRODUCT
inline real dot(const Vector& v1, const Vector& v2)
{
    /*!
   * Calculates dot product between two vectors v1 and v2
   */
#ifdef DEBUG
    assert(v1.getLength() == v2.getLength() );
#endif
    register real d(0.0);
    const idx n = v1.getLength();
//#ifdef USES_OPENMP
//#pragma omp for reduction(+:d)
//#endif
    for (idx i = 0 ; i < n ; ++i)
        d+= v1[i] * v2[i];

    return d;
}
inline void axpy(const real a, const Vector& x, Vector& y)
{
    /*!
   * y = y + a*x
   */
#ifdef USES_OPENMP
#pragma omp parallel for
#endif
    for (idx i = 0 ; i < x.getLength() ; ++i)
    {
        y[i]+= a*x[i];
    }
}// end void axpy

inline void aypx(const real a, Vector& y , const Vector& x)
{
    /*!
   * y = a*y + x
   */
    const idx n = y.getLength();
#ifdef USES_OPENMP
#pragma omp parallel for
#endif
    for (idx i = 0; i < n ; ++i)
    {
        y[i] = a*y[i] + x[i];
    }

}// end void aypx


// VECTOR SCALAR PRODUCT v2 = a*v1
//void multiply(const Vector& v1, const real a, Vector& v2);


// ESTIMATES ERROR NORM FOR SYSTEM Ax=b
// RETURNING (b - Ax).(b - Ax)
real errorNorm2(const IRCMatrix& A, const Vector& x, const Vector& b);
real errorNorm2(const TDMatrix& A, const Vector& x, const Vector& b);

inline real norm(const Vector& x)
{
    return x.getNorm();
}

} // end namespace SpaMtrix

#endif

