#ifndef IRCMATRIX
#define IRCMATRIX
#include <iostream>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <cstring>
#include <omp.h>

#include <setup.h>
#include <vector.h>
#include <fleximatrix.h>

// Interleaved Row Compressed Matrix
class IRCMatrix
{
  /*!
   * Sparse matrix class, using Interleaved Row Compressed storage scheme
   */
protected:
    idx* rows;      // ROW COUNTER
    IndVal* cvPairs;// COLUMN-VALUE PAIRS
    idx nnz;        // NUMBER OF NON-ZEROS
    idx numRows, numCols;
    
    idx getIndex(const idx row, const idx col) const;
public:
    
    IRCMatrix():
        rows(0), cvPairs(0), nnz(0), numRows(0), numCols(0) {}
    IRCMatrix(  const idx numRows, const idx numCols,
                const idx nnz, 
                idx * const rows, IndVal *const cvPairs):
                rows(rows), cvPairs(cvPairs), nnz(nnz), 
                numRows(numRows) , numCols(numCols){}

    IRCMatrix(const IRCMatrix &m);
    IRCMatrix& operator=(const IRCMatrix& m);
    inline IRCMatrix& operator=(const real &s)	
    {
      /*! SETS ALL NONZEROS TO SCALAR s.*/
      for (IndVal *cvp = cvPairs ; cvp < cvPairs + nnz ; cvp++ , cvp->val = s );
      return *this;
    }
    ~IRCMatrix();
    //================================================
    inline idx getnnz()const 
    {
      /*!Returns total number of non-zeros in matrix.*/
      return nnz;
    }
    inline idx getNumRows()const 
    { /*! Returns number of rows in matrx.*/
      return numRows;
    }
    inline idx getNumCols()const 
    {/*!Returns number of columns in matrix*/
      return numCols;
    }
    void sparse_set(const idx row, const idx col , const real val );
    void sparse_add(const idx row, const idx col , const real val );
    real sparse_get(const idx row, const idx col ) const;
    
    // RETURNS VALUE AT (ROW,COL), EVEN IF IT IS ZERO
    real getValue(const idx row, const idx col)const;

    // MATHS OPERATORS
    Vector operator*(const Vector& x) const; // MATRIX VECTOR MULTIPLICATION

    // RETURNS TRUE IS ROW,COL ENTRY EXISTS AS A NON-ZERO. ACTUAL VALUE IS RETURNED IN val
    bool isNonZero(const idx row, const idx col) const;
    bool isNonZero(const idx row, const idx col, real& val) const;
    //===============================================
    // FRIEND FUNCTIONS THAT REQUIRE ACCESS TO PRIVATE
    // DATA FOR PERFORMANCE - THIS MAY GET MESSY
    //===============================================
    // SPAMTRIX-BLAS
    friend void multiply(const IRCMatrix& A, const Vector& x, Vector& b);
    friend real multiply_dot(const IRCMatrix& A, const Vector& x, Vector& b);
    // RELAXERS
    friend void jacobi(const IRCMatrix& A, Vector &x, const Vector &b,const idx &maxIter);
    friend void SOR(const IRCMatrix& A, Vector &x, const Vector& b, idx maxIter);
    
    friend class FlexiMatrix;
    //================================================
    // DEBUG FUNCTIONS

    void spy()const;
    void print() const;

    
    
    
    
};

// MAKES SIMPLE IRCMatrix OUT OF AN ARRAY OF CONNECTED INDEXES
//IRCMatrix& makeIRCMatrix(const idx* t, const idx npt, const idx nt);

#endif
