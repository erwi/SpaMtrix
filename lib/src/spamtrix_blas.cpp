#include <spamtrix_setup.hpp>
#include <spamtrix_blas.hpp>
#include <spamtrix_ircmatrix.hpp>
#include <spamtrix_vector.hpp>
#include <spamtrix_tdmatrix.hpp>
namespace SpaMtrix {


// MATRIX VECTOR MULTIPLICATION b = Ax
void multiply(const IRCMatrix &A,
              const Vector &x,
              Vector &b) {
    /*!
    * Matrix-Vector multiplication Ax=b.
    * \n
    * In python wrappers this function is renamed to IRCMatMul.
    */
#ifdef DEBUG
    assert(A.getNumCols() == x.getLength());
    assert(b.getLength() == x.getLength());
#endif
    // FOR EACH ROW
    const idx n = A.getNumRows();
#ifdef USES_OPENMP
    #pragma omp parallel for
#endif
    for (idx i = 0; i < n; ++i) {
        // FOR EACH COLUMN
        real r(0);
        const idx row_start = A.rows[i];
        const idx row_end   = A.rows[i + 1];
        for (idx j = row_start; j < row_end; ++j) {
            const idx col = A.cvPairs[j].ind;
            r += A.cvPairs[j].val * x[col];
        }//j
        b[i] = r;
    }// i
}

void multiply(const TDMatrix &A,
              const Vector &x,
              Vector &b) {
    /*!
    * Matrix vector muliplication Ax = b
    * where A is a tridiagonal matrix \n \n
    * In python this function is renamed to TDMatMul.
    */
#ifdef DEBUG
    assert(A.size == x.getLength());
#endif
    idx n = b.getLength();
    for (idx i = 0; i < n; i++) {
        b[i] = A.diagonal[i] * x[i];
        if (i > 0)
            b[i]   += A.lower[i - 1] * x[i - 1];
        if (i < (n - 1))
            b[i] += A.upper[i] * x[i + 1];
    }// end for i
}
//
//==============================================================
//
real multiply_dot(const IRCMatrix &A,
                  const Vector &x,
                  Vector &b) {
    /*!
    * MATRIX VECTOR MULTIPLICATION Ax=b
     also returns dot product of x an b. Used in IterativeSolvers
    */
#ifdef DEBUG
    assert(A.getNumCols() == x.getLength());
    assert(b.getLength() == x.getLength());
#endif
    // FOR EACH ROW
    const idx n = A.getNumRows();
    real dp(0);
#ifdef USES_OPENMP
    #pragma omp parallel for reduction(+:dp)
#endif
    for (idx i = 0; i < n; ++i) {
        // FOR EACH COLUMN
        real r(0);
        const idx row_start = A.rows[i];
        const idx row_end   = A.rows[i + 1];
        for (idx j = row_start; j < row_end; ++j) {
            const idx col = A.cvPairs[j].ind;
            r += A.cvPairs[j].val * x[col];
        }//j
        b[i] = r;
        dp += b[i] * x[i];
    }//
    return dp;
}

void scale(const real a, Vector &v) {
    /*!v*= a*/
    idx len = v.getLength();
    for (idx i = 0; i < len; ++i)
        v[i] *= a;
}

real dot(const Vector &v1, const Vector &v2) {
    /*!
    * Calculates dot product between two vectors v1 and v2
    */
#ifdef DEBUG
    assert(v1.getLength() == v2.getLength());
#endif
    real d(0.0);
    const idx n = v1.getLength();
    for (idx i = 0; i < n; ++i)
        d += v1[i] * v2[i];
    return d;
}

void axpy(const real a, const Vector &x, Vector &y) {
    /*!y = y + a*x*/
    for (idx i = 0; i < x.getLength(); ++i)
        y[i] += a * x[i];
}

void aypx(const real a, Vector &y , const Vector &x) {
    /*!y = a*y + x */
    const idx n = y.getLength();
#ifdef USES_OPENMP
    #pragma omp parallel for
#endif
    for (idx i = 0; i < n; ++i)
        y[i] = a * y[i] + x[i];
}

real norm(const Vector &x) {
    return x.getNorm();
}

real errorNorm2(const IRCMatrix &A,
                const Vector &x,
                const Vector &b) {
    /*!
     * Calculates error nor for system Ax = b, using
     * error = (Ax-b).(Ax-b)
     */
    real error(-1);
    Vector temp(x.getLength());
    multiply(A, x, temp);     // temp = Ax
    temp -= b;
    error = dot(temp, temp);
    return error;
}

real errorNorm2(const TDMatrix &A,
                const Vector &x,
                const Vector &b) {
    /*!
     * Calculates error nor for system Ax = b, using
     * error = (Ax-b).(Ax-b)
     */
    real error(-1);
    Vector temp(x.getLength());
    multiply(A, x, temp);     // temp = Ax
    temp -= b;
    error = dot(temp, temp);
    return error;
}

Vector abs(const SpaMtrix::Vector &vin) {
    /*! Returns a copy of vin where all component values are positive*/
    SpaMtrix::Vector vout(vin);
    for (idx i = 0; i < vout.getLength(); ++i)
        vout(i) = fabs(vout(i));
    return vout;
}



} // end namespace SpaMtrix
