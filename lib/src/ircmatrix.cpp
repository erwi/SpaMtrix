#include <iostream>
#include <stdio.h>
#include <assert.h>
#include <algorithm>
#include <stdlib.h>
#include <cstring>
#include <spamtrix_ircmatrix.hpp>
#include <spamtrix_fleximatrix.hpp>
#include <spamtrix_vector.hpp>
namespace SpaMtrix {
IRCMatrix::IRCMatrix(): rows(0), cvPairs(0), nnz(0), numRows(0), numCols(0) {
    /*!Constructs an empty matrix*/
}
IRCMatrix::IRCMatrix(const idx numRows, const idx numCols,
                     const idx nnz,
                     idx *const rows, IndVal *const cvPairs):
    rows(rows), cvPairs(cvPairs), nnz(nnz),
    numRows(numRows) , numCols(numCols) {
    /*!Constructs a sparse matrix where sparsity pattern is
    defined in the arrays 'rows' and 'cvPairs'. The matrix takes ownership of these arrays and releases
    the memory when destructed (i.e. does not make copies).
    */
}
IRCMatrix::IRCMatrix(const IRCMatrix &m):
    rows(NULL),
    cvPairs(NULL),
    nnz(0),
    numRows(0), numCols(0) {
    nnz = m.nnz;
    numRows = m.numRows;
    numCols = m.numCols;
    if (numRows > 0) { // IF NOT EMPTY MATRIX
        rows = new idx [numRows + 1];
        if (!rows) {
            std::cerr << "error in " << __func__ << " could not allocate rows" << std::endl;
            exit(1);
        }
        cvPairs = new IndVal[nnz];
        if (!cvPairs) {
            std::cerr << "error in " << __func__ << " could not allocate cvPairs" << std::endl;
            exit(1);
        }
        memcpy(rows, m.rows, (numRows + 1) * sizeof(idx));
        memcpy(cvPairs, m.cvPairs, nnz * sizeof(IndVal));
    }
}

IRCMatrix::IRCMatrix(IRCMatrix &&m) {
    /*!
     * Copy constructor from rvalue refrence
     */
    // copy data references from m to this
    this->rows = m.rows;
    this->cvPairs = m.cvPairs;
    this->nnz = m.nnz;
    this->numCols = m.numCols;
    this->numRows = m.numRows;
    // Clear m data
    m.rows = nullptr;
    m.cvPairs = nullptr;
}


IRCMatrix::IRCMatrix(const FlexiMatrix &M): rows(NULL), cvPairs(NULL),
    nnz(0), numRows(0), numCols(0) {
    copyFrom(M);
}
IRCMatrix &IRCMatrix::operator=(const IRCMatrix &M) {
    /*!
      Asigment operator. Sets this matrix equal to matrix M. Reallocates memory if necessary
    */
    if (&M == this)
        return *this;
    numRows = M.numRows;
    numCols = M.numCols;
    // IF OTHER MATRIX SIZE IS DIFFERENT FROM CURRENT, MUST REALLOCATE
    if (nnz != M.nnz) {
        nnz = M.nnz;
        // TODO: ADD CHECKS FOR ALLOCATION FAILS
        delete [] rows;
        if (numRows > 0) { // CHECK FOR EMPTY MATRIX
            rows = new idx[numRows + 1];
            memcpy(rows, M.rows, (numRows + 1)*sizeof(idx));
        } else {
            rows = NULL;
        }
        delete [] cvPairs;
        if (nnz > 0) { // CHECK FOR EMPTY MATRIX
            cvPairs = new IndVal[nnz];
            memcpy(cvPairs, M.cvPairs, nnz * sizeof(IndVal));
        } else {
            cvPairs = NULL;
        }
    }
    return *this;
}

IRCMatrix &IRCMatrix::operator=(const real &s) {
    /*! SETS ALL NONZEROS TO SCALAR s.*/
    for (idx i = 0 ; i < nnz ; ++i) {
        cvPairs[i].val = s;
    }
    return *this;
}
IRCMatrix &IRCMatrix::operator=(const FlexiMatrix &M) {
    copyFrom(M);
    return *this;
}
const IRCMatrix IRCMatrix::operator*(const real &s) const {
    /*! Matrix scalar multiplication. Returns a scaled version of self.*/
    // CREATE NEW SPARSE MATRIX OF SAME SIZE AS SELF
    idx *rows_n = new idx[numRows + 1];  // NEW ROW INDEXES
    IndVal *cvPairs_n = new IndVal[nnz]; // NEW COLUMNS/VALUES
    memcpy(rows_n, rows, (numRows + 1)*sizeof(idx));
    for (idx i = 0 ; i < nnz ; ++i) { // FILL NEW COL/VALS
        cvPairs_n[i] = cvPairs[i];
        cvPairs_n[i].val *= s;
    }
    return IRCMatrix(numRows, numCols, nnz, rows_n, cvPairs_n);
}

void IRCMatrix::add(const IRCMatrix &other, const real &scalar) {
#ifdef DEBUG
    assert(numRows >= other.numRows);
    assert(numCols >= other.numCols);
#endif

    // for each row in other matrix
    for (idx i = 0 ; i < other.numRows ; i++) {
        // for each column in other matrix
        const idx row_start = other.rows[i];
        const idx row_end   = other.rows[i + 1];
        if (row_start == row_end) { // this row is empty, i.e. no non-zeros exist on this row
          continue;
        }

        for (idx j = row_start ; j < row_end ; j++) {
            const idx col = other.cvPairs[j].ind;
            const real val = other.cvPairs[j].val;
            sparse_add(i, col, scalar * val);
        }
    }
}

IndVal & IRCMatrix::find(const idx row, const idx col) {
  idx rowStart = this->rows[row];
  idx rowEnd = this->rows[row + 1];

  // binary search between rowStart and rowEnd to find cvPair with column index col
  IndVal *itr = std::lower_bound(&this->cvPairs[rowStart], &this->cvPairs[rowEnd], col,
                                 [](const IndVal &iv, const idx &col) { return iv.ind < col; });

  // return a ref to the found element if it exists
  if (itr != &this->cvPairs[rowEnd] && itr->ind == col) {
    return *itr;
  } else {
    throw std::runtime_error("Index not found (row, col)=(" + std::to_string(row) + ", " + std::to_string(col) + ")");
  }
}

const IndVal& IRCMatrix::find(const idx row, const idx col) const {
  IndVal *begin = &cvPairs[rows[row]]; // pointer to first in row
  IndVal *end = &cvPairs[rows[row + 1]]; // pointer to first in row+1

  // binary search between rowStart and rowEnd to find cvPair with column index col
  IndVal *itr = std::lower_bound(begin, end, col,
                                 [](const IndVal &iv, const idx &col) { return iv.ind < col; });

  // return a ref to the found element if it exists
  if (itr->ind == col && itr != end) {
    return *itr;
  } else {
    throw std::runtime_error("Index not found (row, col)=(" + std::to_string(row) + ", " + std::to_string(col) + ")");
  }
}

IRCMatrix::~IRCMatrix() {
    clear();
}

void IRCMatrix::clear() {
    /*!
     * Clears all data and deallocates memory.
     */
    delete [] rows;
    delete [] cvPairs;
    rows = NULL;
    cvPairs = NULL;
    numCols = 0;
    numRows = 0;
    nnz = 0;
}

void IRCMatrix::copyFrom(const FlexiMatrix &A) {
  if (nnz > 0) {
    clear();
  }
  nnz = A.calcNumNonZeros();
  numRows = A.getNumRows();
  numCols = A.getNumCols();
  cvPairs = new IndVal[nnz];
  if (numRows > 0) { // IF NOT EMPTY MATRIX
    rows = new idx[numRows + 1];
    idx numNonZero = 0;

    // copy data from A to this
    for (idx r = 0; r < numRows; ++r) {
      rows[r] = numNonZero; // index to start of row
      for (auto nz : A.row(r)) { // for each non-zero in row r
        cvPairs[numNonZero] = nz;
        numNonZero++;
      }
    }
    rows[numRows] = numNonZero;
  }
}

void IRCMatrix::sparse_set(const idx row, const idx col, const real val) {
  IndVal &iv = find(row, col);
#ifdef USES_OPENMP
    #pragma omp critical
#endif
  iv.val = val;
}

void IRCMatrix::sparse_add(const idx row, const idx col, const real val) {
  IndVal &iv = find(row, col);
#ifdef USES_OPENMP
    #pragma omp atomic
#endif
  iv.val += val;
}

real IRCMatrix::sparse_get(const idx row, const idx col) const {
  return find(row, col).val;
}

real IRCMatrix::getValue(const idx row, const idx col) const {
    /*!
      Returns value at (row,col). If a non-zero does not exist at
      (row,col), a zero is returned. \n\n
      Note:\n
      Use bool IRCMatrix::isNonZero(row, col, val) instead if it
      is important to know wheteher position (row,col) is a non-zero
      */
    real val;
    return isNonZero(row, col, val) ? val : 0;
}

real* IRCMatrix::getValuePtr(const idx row, const idx col) {
#ifndef NDEBUG
    assert(row < this->numRows);
    assert(col < this->numCols);
#endif

    IndVal *begin = &cvPairs[rows[row]]; // pointer to first in row
    IndVal *end = &cvPairs[rows[row + 1]]; // pointer to first in row+1
    IndVal *itr = std::lower_bound(begin, end, col,
                           [](const IndVal & iv1, const idx& colind) { return iv1.ind < colind; });

    if (itr->ind == col && itr != end) {
        return &itr->val;
    } else {
        return nullptr;
    }
}

  real* IRCMatrix::getValuePtr(const idx row, const idx col) const {
#ifndef NDEBUG
    assert(row < this->numRows);
    assert(col < this->numCols);
#endif

    IndVal *begin = &cvPairs[rows[row]]; // pointer to first in row
    IndVal *end = &cvPairs[rows[row + 1]]; // pointer to first in row+1
    IndVal *itr = std::lower_bound(begin, end, col,
                                   [](const IndVal & iv1, const idx& colind) { return iv1.ind < colind; });

    if (itr->ind == col && itr != end) {
      return &itr->val;
    } else {
      return nullptr;
    }
  }

bool IRCMatrix::isNonZero(const idx row, const idx col) const {
#ifdef DEBUG
    assert(row < this->numRows);
    assert(col < this->numCols_);
#endif
    return getValuePtr(row, col) != nullptr;
}

bool IRCMatrix::isNonZero(const idx row, const idx col, real &val) const {
#ifdef DEBUG
    assert(row < this->numRows);
    assert(col < this->numCols_);
#endif

  auto ptr = getValuePtr(row, col);
  if (ptr != nullptr) {
    val = *ptr;
    return true;
  } else {
    return false;
  }
}

void IRCMatrix::operator*=(const real &s) {
    /*! Modifies matrix by multiplying all values by scalar coefficient s.*/
    for (idx i = 0 ; i < nnz ; i++)
        cvPairs[i].val *= s;
}

Vector IRCMatrix::operator *(const Vector &x) const {
    /*!
     * MATRIX VECTOR MULTIPLICATION Ax=b.
     * reference to b is returned
     */
#ifdef DEBUG
    if (numCols_ != x.getLength()) {
        std::cerr << "error in " << __func__ << " column count is :" << numCols_ <<
                  ", vector length is :" << x.getLength() << std::endl;
        exit(1);
    }
#endif
    Vector b(x.getLength());
    // FOR EACH ROW
    for (idx i = 0 ; i < getNumRows() ; i++) {
        // FOR EACH COLUMN
        real r(0);
        const idx row_start = rows[i];
        const idx row_end   = rows[i + 1];
        for (idx j = row_start ; j < row_end ; j++) {
            const idx col = cvPairs[j].ind;
            r += cvPairs[j].val * x[col];
        }// end for jj
        b[i] = r;
    }//
    return b;
}


// #ifdef DEBUG


void IRCMatrix::spy() const {
    std::cout << std::endl;
    std::cout << "IRCMatrix size = " << this->numRows << " , " << this->numCols << " nnz = " << this->nnz << std::endl;
    for (idx row = 0 ; row < this->numRows ; row++) {
        for (idx col = 0 ; col < this->numCols ; col++) {
            char marker;
            if (this->isNonZero(row, col))
                marker = '#';
            else
                marker = '.';
            printf("%c", marker);
        }
        printf("\n");
    }
}

void IRCMatrix::print() const {
    std::cout << std::endl;
    std::cout << "IRCMatrix size = " << this->numRows << " , " << this->numCols << " nnz = " << this->nnz << std::endl;
    for (idx row = 0 ; row < this->numRows ; row++) {
        for (idx col = 0 ; col < this->numCols ; col++) {
            real val = this->getValue(row, col);
            printf("%1.6f\t", val);
        }
        printf("\n");
    }
}

} // end namespace SpaMtrix
// #endif
