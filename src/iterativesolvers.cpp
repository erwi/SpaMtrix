#include <iterativesolvers.h>
#include <spamtrix_blas.h>
#include <densematrix.h>




namespace SpaMtrix
{

IterativeSolvers::IterativeSolvers():
		  maxIter(0),
		  maxInnerIter(0),
		  toler(0)
{}

IterativeSolvers::IterativeSolvers(const idx maxIter,
				   const real toler):
				   maxIter(maxIter),
				   maxInnerIter(0),
				   toler(toler)
{}
IterativeSolvers::IterativeSolvers(const idx maxIter,
				   const idx maxInnerIter,
				   const real toler):
				   maxIter(maxIter),
				   maxInnerIter(maxInnerIter),
				   toler(toler)
{}

bool IterativeSolvers::pcg(const IRCMatrix &A,
                           Vector &x,
                           const Vector &b,
                           const Preconditioner &M)
{
    /*!
      Solves Ax=b using the preconditioned conjugate gradient method.
      */
    const idx N = x.getLength();
    real resid(100.0);
    Vector p(N), z(N), q(N);

    real alpha;
    real normr(0);
    real normb = norm(b);
    real rho(0), rho_1(0), beta(0);

    Vector r = b - A*x;

    if (normb == 0.0)
        normb = 1;
    resid = norm(r) / normb;
    if (resid <= IterativeSolvers::toler) 
    {
	IterativeSolvers::toler = resid;
        IterativeSolvers::maxIter = 0;
        return true;
    }
    // MAIN LOOP
    idx i = 1;
    for (; i <= IterativeSolvers::maxIter; i++)
    {
  
        M.solveMxb(z,r);
        rho = dot(r, z);

        if (i == 1)
            p = z;
        else
        {
            beta = rho / rho_1;
            aypx(beta, p, z); // p = beta*p + z;
        }

        // CALCULATES q = A*p AND dp = dot(q,p)
        real dp = multiply_dot(A,p,q);
        alpha = rho / dp;
        normr = 0;
#ifdef USES_OPENMP	
#pragma omp parallel for reduction(+:normr)
#endif
        for (idx j = 0 ; j < N ; ++j)
        {
            x[j] += alpha*p[j];// x + alpha(0) * p;
            r[j] -= alpha*q[j];// r - alpha(0) * q;
            normr += r[j]*r[j];
        }
        normr = sqrt(normr);
        resid = normr/normb;
        if (resid <= IterativeSolvers::toler)
        {
            IterativeSolvers::toler = resid;
            IterativeSolvers::maxIter = i;
	    return true;
        }
        rho_1 = rho;
    }
    IterativeSolvers::toler = resid;
    return false;
}

//=====================================================
// FUNCTIONS USED BY GMRES ONLY
inline void GeneratePlaneRotation(const real &dx, const real &dy, real &cs, real &sn)
{
    if (dy == 0.0)
    {
        cs = 1.0;
        sn = 0.0;
    } else if (fabs(dy) > fabs(dx))
    {
        real temp = dx / dy;
        sn = 1.0 / sqrt( 1.0 + temp*temp );
        cs = temp * sn;
    } else
    {
        real temp = dy / dx;
        cs = 1.0 / sqrt( 1.0 + temp*temp );
        sn = temp * cs;
    }
}



inline void ApplyPlaneRotation(real &dx, real &dy, const real &cs, const real &sn)
{
    const real temp  =  cs * dx + sn * dy;
    dy = -sn * dx + cs * dy;
    dx = temp;
}


inline void Update(Vector &x,
                   const idx k,
                   const DenseMatrix &h,
                   Vector &temp,
                   const Vector v[]
                   )
{
    // Backsolve:
    for (int i = k; i >= 0; --i)
    {
        temp(i) /= h(i,i);
        for (int j = i - 1; j >= 0; --j)
            temp(j) -= h(j,i) * temp(i);
    }

    for (idx j = 0; j <= k; ++j)
        x += v[j] * temp(j);
}


//template < class real >
inline real abs(const real x)
{
    return (x > 0 ? x : -x);
}


bool IterativeSolvers::gmres(const IRCMatrix &A,
                             Vector &x,
                             const Vector &b,
                             const Preconditioner &M)
{
    const idx N = x.getLength();
    idx i, j = 1, k;
    Vector s(maxInnerIter+1);
    Vector cs(maxInnerIter+1);
    Vector sn(maxInnerIter+1);
    Vector w(N);

    real normb = norm(M.solve(b));
    Vector r = M.solve(b - A * x);
    real beta = norm(r);

    if (normb == 0.0)
        normb = 1;
    real res(norm(r) / normb);
    if (res <= toler)
    {
        toler = res;
        maxIter = 0;
        return true;
    }

    Vector *v = new Vector[maxInnerIter+1];
    for (idx id = 0; id < maxInnerIter+1; ++id)
        v[id] = Vector(N);


    // CREATE HESSENBERG MATRIX NEEDED TO STORE INTERMEDIATES
    DenseMatrix H(maxInnerIter+1, maxInnerIter);

    Vector temp(N);
    Vector temp2(maxInnerIter+1);

    // MAIN LOOP
    while (j <= maxIter)
    {
        v[0] = r * (1.0 / beta);
        s = 0.0;

        s(0) = beta;

        //cout << "j :" << j << endl;
        // INNER ITERATIONS
        for (i = 0; i < maxInnerIter && j <= maxIter; i++, j++)
        {

            // CALCULATE w = M^{-1}(A*v[i])
            multiply(A,v[i], temp);
            M.solveMxb(w,temp);

            // PRE-CALCULATE DOT PRODUCTS IN PARALLEL
            // H(k,i) = dot( v[k], w)
#ifdef USES_OPENMP
            #pragma omp parallel for
#endif
            for (k = 0; k <= i ; ++k)
            {
                register real dp(0);
                for (idx id = 0 ; id < N ; ++id)
                    dp+= w[id]*v[k][id];

                H(k,i) = dp;//dot(w,v[k]);
            }


            for (k = 0; k <= i; ++k)
            {   // w -= v[k]*H(k,i) without temporaries
                register real tempr = H(k,i);
#ifdef USES_OPENMP
                #pragma omp parallel for // why is this loop so critical??
#endif
                for (idx id = 0 ; id < N ; ++id)
                    w[id] -= v[k][id]*tempr;
            }

            // BELOW PARALLEL REGION CALCULATES:
            // H(i+1,i) = norm(w);
            // v[i+1] = w * (1.0 / H(i+1, i));
            H(i+1,i) = 0;
            real tempr(0);
#ifdef USES_OPENMP
            #pragma omp parallel shared(tempr)
#endif
            {
#ifdef USES_OPENMP
                #pragma omp for reduction(+:tempr)
#endif
                for (idx id = 0 ; id < N ; ++id)
                    tempr += w[id]*w[id]; //norm(w);
#ifdef USES_OPENMP
                #pragma omp single
#endif
                {
                    H(i+1,i) = sqrt(tempr );
                    tempr = (1.0/H(i+1,i) );
                }
#ifdef USES_OPENMP
                #pragma omp for
#endif
                for (idx id = 0 ; id < N ; ++id)
                    v[i+1][id] = w[id]*tempr;
            }// end for omp parallel



            for (k = 0; k < i; k++)
                ApplyPlaneRotation(H(k,i), H(k+1,i), cs(k), sn(k));

            GeneratePlaneRotation(H(i,i), H(i+1,i), cs(i), sn(i));
            ApplyPlaneRotation(H(i,i), H(i+1,i), cs(i), sn(i));
            ApplyPlaneRotation(s(i), s(i+1), cs(i), sn(i));

            res = fabs(s(i+1) ) / normb;

            if (res < toler)
            {

                // COPY S INTO temp WITHOUT RESIZING
                for (idx id = 0 ; id < maxInnerIter+1 ; ++id)
                    temp2[id] = s[id];
                Update(x, i, H, temp2, v);
                toler = res;
                maxIter = j;
                delete [] v;
                return true;
            }

        }// end for i IINNER ITERATIONS

        // COPY S INTO temp WITHOUT RESIZING
        for (idx id = 0 ; id < maxInnerIter+1 ; ++id)
            temp2[id] = s[id];

        Update(x, maxInnerIter - 1, H, temp2, v);

        //multiply(A, x, temp);     //r = M.solve(b - A * x);
        M.solveMxb(r, b-A*x);


        beta = norm(r);
        res = beta / normb;
        if (res < toler)
        {
            toler = res;
            maxIter = j;
            delete [] v;
            return true;
        }
    }

    toler = res;
    delete [] v;
    return false;
}
}// end namespace SpaMtrix