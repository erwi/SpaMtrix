

#ifndef RCMATRIX_H
#define RCMATRIX_H

#include "mvvd.h" // DOUBLE VECTOR
#include "mvvi.h" // INTEGER VECTOR

//#include "vecdefs.h"
//#include VECTOR_H

//class CompCol_Mat_double;
//class Coord_Mat_double;

class RCMatrix {

private:
       Vector val_;       // data values (nz_ elements)
       MV_Vector_int    rowptr_;    // row_ptr (dim_[0]+1 elements)
       MV_Vector_int    colind_;    // col_ind  (nz_ elements)

       int base_;                 // index base: offset of first element
       int nz_;                   // number of nonzeros
       int dim_[2];               // number of rows, cols
  
public:
       RCMatrix(void);
       RCMatrix(const RCMatrix &S);

      // RCMatrix(const Coord_Mat_double &CO);
       RCMatrix(int M, int N, int nz, double *val, int *r,
                               int *c, int base=0);
       RCMatrix(int M, int N, int nz,
                              const Vector &val,
                              const MV_Vector_int &r, 
			      const MV_Vector_int &c,
                              int base=0);
      ~RCMatrix() {};
    
/*******************************/
/*  Access and info functions  */
/*******************************/

       double&      val(int i) { return val_(i); }
       int&         row_ptr(int i) { return rowptr_(i); }
       int&         col_ind(int i) { return colind_(i);}

       const double&      val(int i) const { return val_(i); }
       const int&         row_ptr(int i) const { return rowptr_(i); }
       const int&         col_ind(int i) const { return colind_(i);}

       int dim(int i) const {return dim_[i];}
       int size(int i) const {return dim_[i];}
       int NumNonzeros() const {return nz_;}
       int base() const {return base_;}

/*******************************/
/*  Assignment  operator       */
/*******************************/

       RCMatrix& operator=(const RCMatrix &R);
       RCMatrix& newsize(int M, int N, int nz);

/***********************************/
/*  General access function (slow) */
/***********************************/

       double       operator() (int i, int j) const;        
       double&      set(int i, int j);

/***********************************/
/*  Matrix/Vector multiply         */
/***********************************/

       Vector operator*(const Vector &x) const;
       Vector trans_mult(const Vector &x) const;

};

//std::ostream& operator << (std::ostream & os, const RCMatrix & mat);

#endif  
/* RCMatrix_H */

