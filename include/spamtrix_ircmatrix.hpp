#ifndef IRCMATRIX
#define IRCMATRIX

#include <spamtrix_setup.hpp>
namespace SpaMtrix
{
// forward decarations
class FlexiMatrix;
class Vector;

// Interleaved Row Compressed Matrix
class IRCMatrix{
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

    IRCMatrix();
    IRCMatrix(  const idx numRows, const idx numCols,
                const idx nnz,
                idx * const rows, IndVal *const cvPairs);
    IRCMatrix(const IRCMatrix &m);
    IRCMatrix(const FlexiMatrix &M);
    IRCMatrix& operator=(const IRCMatrix& m);
    IRCMatrix& operator=(const real &s);
    IRCMatrix& operator=(const FlexiMatrix &m);
    virtual ~IRCMatrix();
    //================================================
    void clear();
    void copyFrom(const FlexiMatrix& A); // reallocates using data from fleximatrix
    idx getnnz()const;      // RETURNS NUMBER OF NONZEROS
    idx getNumRows()const;  // RETURNS MATRIX ROW COUNT
    idx getNumCols() const; // RETURNS MATRIX COLUMN COUNT

    void sparse_set(const idx row, const idx col , const real val );
    void sparse_add(const idx row, const idx col , const real val );
    real sparse_get(const idx row, const idx col ) const;

    // RETURNS VALUE AT (ROW,COL), EVEN IF IT IS ZERO
    real getValue(const idx row, const idx col)const;

    // MATHS OPERATORS
    Vector operator*(const Vector& x) const;    // MATRIX VECTOR MULTIPLICATION
    void operator*=(const real &s);             // IN-PLACE MULTIPLICATION BY A SCALAR COEFFICIENT
    const IRCMatrix operator*(const real &s)const;  // RETURNS A SCALED VERSION OF THE MATRIX
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
    friend class FlexiMatrix;
    //================================================
    // DEBUG FUNCTIONS
    void spy()const;
    void print() const;
};

inline const IRCMatrix operator*(const real& s, const IRCMatrix &M)
{
  /*! Implements scalar*Matrix, in terms of Matrix*scalar, to achieve commutativity*/
    return M*s;
}

// IMPLEMENTATIONS OF INLINDED METHODS - NO NEW DECLARATIONS BELOW THIS
inline idx IRCMatrix::getNumCols() const {/*!Returns matrix column count*/ return numCols;}
inline idx IRCMatrix::getNumRows() const {/*!Returns matrix row count*/ return numRows;}
inline idx IRCMatrix::getnnz() const {/*!Returns number of nonzeros in matrix*/ return nnz;}


} // end namespace SpaMtrix
#endif
