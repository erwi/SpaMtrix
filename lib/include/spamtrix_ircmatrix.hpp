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

    IndVal & find(const idx row, const idx col);
    [[nodiscard]] const IndVal& find(const idx row, const idx col) const;
public:

    IRCMatrix();
    IRCMatrix(  const idx numRows, const idx numCols,
                const idx nnz,
                idx * const rows, IndVal *const cvPairs);
    IRCMatrix(const IRCMatrix &m);
    IRCMatrix(IRCMatrix &&m);
    IRCMatrix(const FlexiMatrix &M); // copy constructor with move semantics
    IRCMatrix& operator=(const IRCMatrix& m);
    IRCMatrix& operator=(const real &s);
    IRCMatrix& operator=(const FlexiMatrix &m);
    virtual ~IRCMatrix();
    //================================================
    void clear();
    /** Makes a copy of FlexiMatrix A */
    void copyFrom(const FlexiMatrix& A);
    idx getnnz()const;      // RETURNS NUMBER OF NONZEROS
    idx getNumRows()const;  // RETURNS MATRIX ROW COUNT
    idx getNumCols() const; // RETURNS MATRIX COLUMN COUNT

    void sparse_set(const idx row, const idx col , const real val );
    void sparse_add(const idx row, const idx col , const real val );
    real sparse_get(const idx row, const idx col ) const;

    /** Returns value at (row, col), even if it is zero */
    real getValue(const idx row, const idx col) const;

    /** Returns a pointer to the value at (row, col)if it exists, or nullptr otherwise */
    [[nodiscard]] real* getValuePtr(const idx row, const idx col);

    /** Returns a pointer to the value at (row, col)if it exists, or nullptr otherwise */
    [[nodiscard]] real* getValuePtr(const idx row, const idx col) const;

    // MATHS OPERATORS
    /*! MATRIX VECTOR MULTIPLICATION Ax=b.
    * reference to b is returned
    */
    Vector operator*(const Vector& x) const;    // MATRIX VECTOR MULTIPLICATION
    void operator*=(const real &s);             // IN-PLACE MULTIPLICATION BY A SCALAR COEFFICIENT
    const IRCMatrix operator*(const real &s) const;  // RETURNS A SCALED VERSION OF THE MATRIX

    /**
     * Adds the other matrix, scaled by a scalar, to this matrix. Note that the sparsity patterns must match so that
     * the other matrix can not contain any non-zeros at locations that are zero in this matrix.
     */
    void add(const IRCMatrix& other, const real& scalar = 1.0);

    /** return true if storage exists at given row, col */
    [[nodiscard]] bool isNonZero(const idx row, const idx col) const;
    /** return true if storage exists at given row, col and the value is copied to val */
    [[nodiscard]] bool isNonZero(const idx row, const idx col, real& val) const;
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
