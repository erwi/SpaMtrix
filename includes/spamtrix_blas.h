#ifndef SPAMTRIX_BLAS_H
#define SPAMTRIX_BLAS_H
#include <setup.h>
#include <ircmatrix.h>
#include <vector.h>
#include <assert.h>
#include <iostream>
#include <stdlib.h>
#include <omp.h>
// MATRIX VECTOR MULTIPLICATION b = Ax
void multiply(const IRCMatrix& A,
	     const Vector& x,
	     Vector& b);

// VECTOR DOT PRODUCT
real dot(const Vector& v1, const Vector& v2);

// VECTOR SCALAR PRODUCT v2 = a*v1
void multiply(const Vector& v1, const real a, Vector& v2);


#endif

