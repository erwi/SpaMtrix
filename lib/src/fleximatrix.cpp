#include <spamtrix_fleximatrix.hpp>
#include <algorithm>
#include <assert.h>
#include <iostream>
#include <stdio.h>

namespace SpaMtrix {

FlexiMatrix::FlexiMatrix(const IRCMatrix &A) {
    nonZeros = std::vector<std::vector<IndVal> >(A.getNumRows());
    for (idx r = 0; r < A.getNumRows(); r++) {
        idx rowStart = A.rows[r];
        idx rowEnd   = A.rows[r + 1];
        nonZeros[r] = std::vector<IndVal>(&A.cvPairs[rowStart] , &A.cvPairs[rowEnd]);
    }
}

FlexiMatrix::~FlexiMatrix() { }

idx FlexiMatrix::calcNumNonZeros() const {
    /*!
     * Returns number of non-zeros allocated. Does this by calculating it, so this gets
     * linearly slower with larger matrices
     */
    idx nnz(0);
    for (idx i = 0 ; i < nonZeros.size() ; ++i)
        nnz += nonZeros[i].size();
    return nnz;
}

void FlexiMatrix::addNonZero(size_t row, size_t col, real value) {
    // CHECK DIMENSION1 VECTOR SIZE
#ifdef DEBUG
    if (row >= getNumRows())
        std::cerr << "row = " << row << "num rows = " << getNumRows() << std::endl;
    assert(row < getNumRows());
#endif
    // Append empty rows if needed
    while (getNumRows() <= row) {
        nonZeros.emplace_back();
    }
    // IF EMPTY ROW OR NEW VALUE PLACED AT END, JUST APPEND TO ROW AND EXIT FUNCTION
    if ((nonZeros[row].size() == 0) ||
        (nonZeros[row].back().ind < col)) {
        nonZeros[row].emplace_back(col, value);
        numCols_ = std::max(numCols_, (size_t) col + 1);
        return;
    }
    // FIND CORRECT POSITION BY SEARCHING FOR COLUMN POSITIONS.
    // ALL COLUMN VALUES MUST INCREASE, I.E MAINTAINING AN ASCENDIGLY
    // SORTED VECTOR OF NON-ZERO COLUMNS
    // LOWER BOUND RETURNS ITERATOR TO FIRST ELEMENT THAT DOES NOT
    // COMPARE TO LESS THAN dim2
    IndVal temp(col, value);
    auto itr = std::lower_bound(nonZeros[row].begin(),
                            nonZeros[row].end(),
                            temp,
                          [](const IndVal & iv1, const IndVal & iv2) { return iv1.ind < iv2.ind; }
                        );
    // IF ADDING DUPLICATE NONZERO, ONLY WRITE VALUE TO EXISTING MEMORY
    if (itr->ind == temp.ind) {
        itr->val = temp.val;
    }
    // OTHERWISE ADD NEW VALUE
    else {
        nonZeros[row].insert(itr, temp);
    }

    numCols_ = std::max(numCols_, (size_t) temp.ind + 1);
}

real FlexiMatrix::getValue(size_t row, size_t col) const {
    if (row >= getNumRows()) {
        return 0.;
    }
    IndVal temp(col, 0.0);
    // Find iterator to first element that is not less than col (i.e. is equal or greater)
    auto itr = std::lower_bound(nonZeros[row].begin(), nonZeros[row].end(), temp,
                [](const IndVal & iv1, const IndVal & iv2) { return iv1.ind < iv2.ind; });

    if (itr == nonZeros[row].end() || itr->ind != col) {
        return 0.;
    } else {
        return itr->val;
    }
}

std::vector<IndVal>& FlexiMatrix::row(size_t row) {
  assert(row < nonZeros.size());
  return nonZeros[row];
}

const std::vector<IndVal>& FlexiMatrix::row(size_t row) const {
  assert(row < nonZeros.size());
  return nonZeros[row];
}

void FlexiMatrix::setValue(const idx dim1, const idx dim2, const real val) {
    /*!
      Sets value at (dim1,dim2) to value val. If a non-zero does not already
      exist at (dim1,dim2), a new non-zero is inserted
      */
#ifdef DEBUG
    assert(dim1 < numDim1);
    assert(dim2 < numDim2);
#endif
    real *nnzval;
    if (isNonZero(dim1, dim2, nnzval)) {    // IF NON-ZERO STORAGE ALREADY EXISTS, SET ITS VALUE
        printf("update(%d,%d)=%e\n", dim1, dim2, val);
        *nnzval = val;
    } else                                  // OTHERWISE INSERT NEW NONZERO LOCATION
        addNonZero(dim1, dim2, val);
}


void FlexiMatrix::print() const {
    /*!
    Prints the values held in the matrix to stdout
    */
    printf("FlexiMatrix %zu, %zu\n", getNumRows(), getNumCols());
    for (idx r = 0 ; r < nonZeros.size() ; r++) {
        for (idx c = 0 ; c < getNumCols() ; c++) {
            real val = this->getValue(r, c);
            printf("%1.3f\t", val);
        }
        printf("\n");
    }
    fflush(stdout);
}

bool FlexiMatrix::isNonZero(const idx dim1, const idx dim2, real *&val) {
    /*!
    Returns true if memory is allocated for a non-zero value at index dim1,dim2.

    Input argument val is a pointer to location of addres of stored value if non-zero
    (or NULL pointer otherwise).
    */
    val = NULL; // ASSUME STORAGE IS NOT ALLOCATED AT dim1,dim2
    if (dim1 >= (idx) nonZeros.size())
        return false;
    // IF VECTOR HAS NOT INITIALISED
    if (nonZeros[dim1].empty())
        return false;
    // SET ITERATOR TO FIRST ELEMENT THAT DOES NOT COMPARE TO LESS
    // THAN dim2.
    IndVal iv(dim2, 0.0);
    std::vector<IndVal>::iterator itr =
        std::lower_bound(nonZeros[dim1].begin(),
                         nonZeros[dim1].end(),
                         iv,
    [](const IndVal & iv1, const IndVal & iv2) {
        return iv1.ind < iv2.ind;   // lamda column comparison
    }
                        );
// itr WILL POINT TO END IF ALL INDEXES ARE LESS THAN dim2
    if (itr == nonZeros[dim1].end())
        return false;
// itr NOW POINTS EITHER TO dim2 OR THE NEXT AFTER IT
    if ((*itr).ind == dim2) {
        val = &(itr->val); // SET POINTER TO VALUE AT INDEX dim1,dim2
        return true;
    } else
        return false; // NON-ZERO AT INDEX dim1,dim2 NOT FOUND
}
} // end namespace SpaMtrix
