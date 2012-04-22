#ifndef PCG_H
#define PCG_H
#include <setup.h>
#include <ircmatrix.h>
#include <vector.h>
#include <diagpreconditioner.h>
#include <spamtrix_blas.h>

bool pcg(const IRCMatrix& A,
	 const Vector& b,
	 Vector& x,
	 const DiagPreconditioner& M,
	 idx& maxIter,
	 real& toler);

#endif

