#ifndef SPAMTRIX_BLAS_H
#define SPAMTRIX_BLAS_H


#include <spamtrix_setup.hpp>


namespace SpaMtrix
{
class IRCMatrix;
class TDMatrix;
class Vector;

// MATRIX VECTOR MULTIPLICATION b = Ax
void multiply(const IRCMatrix& A,
              const Vector& x,
              Vector& b);
// TRIADIAGONAL MATRIX MULTIPLICATION
void multiply(const TDMatrix& A,
              const Vector& x,
              Vector& b);
// CONVENIENCE , IN-PLACE SPECIAL MULTIPLICATION...
real multiply_dot(const IRCMatrix& A,
                     const Vector& x,
                     Vector& b);
// RETURNS A POSITIVE COPY OF A VECTOR 
SpaMtrix::Vector abs(const SpaMtrix::Vector &vin);

// ======================================
// BLAS LEVEL 1 FUNCTIONS
// ======================================
void scale(const real a, Vector& v);            // v*= a
real dot(const Vector& v1, const Vector& v2);   // VECTOR DOT PRODUCT
void axpy(const real a, const Vector& x, Vector& y);  // y+= a*x
void aypx(const real a, Vector& y , const Vector& x); // y = a*y + x


// ESTIMATES ERROR NORM FOR SYSTEM Ax=b
// RETURNING (b - Ax).(b - Ax)
real errorNorm2(const IRCMatrix& A, const Vector& x, const Vector& b);
real errorNorm2(const TDMatrix& A, const Vector& x, const Vector& b);
real norm(const Vector& x);
} // end namespace SpaMtrix

#endif

