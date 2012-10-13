
#include <setup.h>
#include <vector.h>
#include <ircmatrix.h>
#include <spamtrix_blas.h>
#include <preconditioner.h>

//template < class Matrix, class Vector, class Real >
bool
cg(const IRCMatrix &A,
   Vector &x,
   const Vector &b,
   const Preconditioner &M,
   idx &max_iter,
   real &tol)
{
    /*!
      Preconditioned conjugate gradient method.
    */


    const idx N = x.getLength();

    Vector p( N );
    Vector z( N );
    Vector q( N );
    Vector r = b - A*x;

    real alpha(0), beta(0), rho(0), rho_1(0), resid(0);
    real normb(0);// = norm(b);
    real normr(0);// = norm(r);

    for (idx j = 0 ; j < N ; ++j)
    {
        normb += b[j]*b[j];
        normr += r[j]*r[j];
    }
    normb = normb == 0 ? 1.0:normb;
    if ((resid = normr/ normb) <= tol) {
        tol = resid;
        max_iter = 0;
        return 0;
    }

    M.solveMxb(z,r); // z = M^{-1} r
    rho = dot(r, z);

    // MAIN LOOP
    for (idx i = 0; i < max_iter; ++i)
    {
        //rho = 0.0;
        for (idx j = 0 ; j < N ; ++j)
        {
           // rho += r[j]*z[j];
            p[j] = z[j]+beta*p[j];
        }

        multiply(A,p,q); // q = A*p

        alpha = rho / dot(p, q);

        // CALCULATE x = x + a*p
        //           r = r - a*q
        //   AND     norm(r), WHICH IS NEEDED LATER
        real normr(0);
        for (idx j = 0 ; j < N ; ++j)
        {
            x[j] += alpha*p[j];
            r[j] -= alpha*q[j];
            normr += r[j]*r[j];
        }


        // IF TOLERANCE ACHIVED, EXIT
        if ((resid = normr / normb) <= tol)
        {
            tol = resid;
            max_iter = i;
            return true;
        }

        // FOR NEXT ITERATION
        rho_1 = rho;
        M.solveMxb(z,r); // z = M^{-1} r
        rho = dot(r, z);
        beta = rho / rho_1;

    }// end for i, MAIN LOOP

    tol = resid;
    return false;
}

