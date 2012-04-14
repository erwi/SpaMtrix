

#include <iostream>
#include <stdlib.h>

#include <rcmatrix.h>
//#include "coord_double.h"

//#include "spblas.h"

/*****************************/
/*  Constructor(s)           */
/*****************************/

RCMatrix::RCMatrix(void)
        : val_(0), rowptr_(0), colind_(0), base_(0), nz_(0)
{
        dim_[0] = 0;
        dim_[1] = 0;
}

/*****************************/
/*  Copy constructor         */
/*****************************/

RCMatrix::RCMatrix(const RCMatrix &S) :
        val_(S.val_), rowptr_(S.rowptr_), colind_(S.colind_), base_(S.base_), 
        nz_(S.nz_)
{
        dim_[0] = S.dim_[0];
        dim_[1] = S.dim_[1];
}

/***********************************/
/*  Construct from storage vectors */
/***********************************/

RCMatrix::RCMatrix(int M, int N, int nz, double *val,
                                     int *r, int *c, int base) :
        val_(val, nz), rowptr_(r, M+1), colind_(c, nz), base_(base), nz_(nz) 
{
        dim_[0] = M;
        dim_[1] = N;
}

RCMatrix::RCMatrix(int M, int N, int nz,
                   const Vector &val, const MV_Vector_int &r,
                   const MV_Vector_int &c, int base) :
        val_(val), rowptr_(r), colind_(c), base_(base), nz_(nz)
{
        dim_[0] = M;
        dim_[1] = N;
} 

/***************************/
/* Assignment operator...  */
/***************************/

RCMatrix& RCMatrix::operator=(const RCMatrix &R)
{
        dim_[0] = R.dim_[0];
        dim_[1] = R.dim_[1];
        base_   = R.base_;
        nz_     = R.nz_;
        val_    = R.val_;
        rowptr_ = R.rowptr_;
        colind_ = R.colind_;
        return *this;
}

/***************************/
/* newsize()               */
/***************************/

RCMatrix& RCMatrix::newsize(int M, int N, int nz)
{
        dim_[0] = M;
        dim_[1] = N;
        nz_     = nz;
        val_.newsize(nz);
        rowptr_.newsize(M+1);
        colind_.newsize(nz);
        return *this;
}

/*********************/
/*   Array access    */
/*********************/


double RCMatrix::operator()(int i, int j)  const
{
        for (int t=rowptr_(i); t<rowptr_(i+1); t++)
           if (colind_(t) == j) return val_(t);
        if (i < dim_[0] && j < dim_[1]) return 0.0;
        else
        {
            std::cerr << "Array accessing exception -- out of bounds." << "\n";
            return (0);
        }
}


double& RCMatrix::set(int i, int j)
{
        for (int t=rowptr_(i); t<rowptr_(i+1); t++)
           if (colind_(t) == j) return val_(t);
        std::cerr << "Array element (" << i << "," << j << 
                                ") not in sparse structure -- cannot assign."
             << "\n";          
        exit(1);
    return val_(0);   // // return to suppress compiler warning message
}


/***************************************/
/*  Matrix-MV_Vector multiplication    */
/***************************************/

Vector RCMatrix::operator*(const Vector &x) const
{
    // NEEDS REIMPLEMENTATION
    /*
    int M = dim_[0];
        int N = dim_[1];

//      Check for compatible dimensions:
        if (x.size() != N) 
        {
           std::cerr << "Error in CompCol Matvec -- incompatible dimensions." 
                << "\n";
           exit(1);
           return x;      // return to suppress compiler warning message
         }

         Vector result(M, 0.0);
         Vector work(M);
  
         int descra[9];
         descra[0] = 0;
         descra[1] = 0;
         descra[2] = 0;

         F77NAME(dcsrmm) (0, M, 1, N, 1.0,
                  descra, &val_(0), &colind_(0), &rowptr_(0),
                  &x(1), N, 1.0, &result(0), M,
                  &work(0), M);

         return result;
         */
}


/*************************************************/
/* Matrix-Transpose-MV_Vector multiplication...  */
/*************************************************/

Vector RCMatrix::trans_mult(const Vector &x)
          const
{
    // NEED REIMPLEMENTATION
    /*
    int M = dim_[0];
          int N = dim_[1];

//        Check for compatible dimensions:
          if (x.size() != M) 
          {
             std::cerr << "Error in CompCol Matvec -- incompatible dimensions." 
                  << "\n";
             exit(1);
             return x;   // return to suppress compiler warning message
           }

           Vector result(N, 0.0);
           Vector work(N);
  
           int descra[9];
           descra[0] = 0;
           descra[1] = 0;
           descra[2] = 0;

           F77NAME(dcsrmm) (1, N, 1, M, 1.0,
                        descra, &val_(0), &colind_(0), &rowptr_(0),
                    &x(0), M, 1.0, &result(1), N,
                    &work(1), N);

           return result;
           */
}


