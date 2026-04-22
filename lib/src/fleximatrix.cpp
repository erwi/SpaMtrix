// Copyright (c) 2026 Eero
// SPDX-License-Identifier: MIT
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
        for (const auto &entry : nonZeros[r]) {
            numCols_ = std::max(numCols_, static_cast<size_t>(entry.ind) + 1);
        }
    }
}

FlexiMatrix::~FlexiMatrix() { }

idx FlexiMatrix::calcNumNonZeros() const {
    idx nnz(0);
    for (idx i = 0 ; i < nonZeros.size() ; ++i)
        nnz += nonZeros[i].size();
    return nnz;
}

void FlexiMatrix::addNonZero(size_t row, size_t col, real value) {
    // Check the row index.
#ifdef DEBUG
    if (row >= getNumRows())
        std::cerr << "row = " << row << "num rows = " << getNumRows() << std::endl;
    assert(row < getNumRows());
#endif
    // Append empty rows if needed.
    while (getNumRows() <= row) {
        nonZeros.emplace_back();
    }
    // Append directly if the row is empty or the new entry belongs at the end.
    if ((nonZeros[row].size() == 0) ||
        (nonZeros[row].back().ind < col)) {
        nonZeros[row].emplace_back(col, value);
        numCols_ = std::max(numCols_, (size_t) col + 1);
        return;
    }
    // Find the correct insertion position while keeping columns sorted.
    IndVal temp(col, value);
    auto itr = std::lower_bound(nonZeros[row].begin(),
                            nonZeros[row].end(),
                            temp,
                          [](const IndVal & iv1, const IndVal & iv2) { return iv1.ind < iv2.ind; }
                        );
    // Update an existing entry if the column already exists.
    if (itr->ind == temp.ind) {
        itr->val = temp.val;
    }
    // Otherwise insert a new entry.
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
    // Find the first entry that is not less than col.
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
#ifdef DEBUG
    assert(dim1 < numDim1);
    assert(dim2 < numDim2);
#endif
    real *nnzval;
    if (isNonZero(dim1, dim2, nnzval)) {
        printf("update(%d,%d)=%e\n", dim1, dim2, val);
        *nnzval = val;
    } else
        addNonZero(dim1, dim2, val);
}


void FlexiMatrix::print() const {
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
    // Assume storage is not allocated at dim1, dim2.
    val = NULL;
    if (dim1 >= (idx) nonZeros.size())
        return false;
    // If the row vector has not been initialized.
    if (nonZeros[dim1].empty())
        return false;
    // Set the iterator to the first element not less than dim2.
    IndVal iv(dim2, 0.0);
    std::vector<IndVal>::iterator itr =
        std::lower_bound(nonZeros[dim1].begin(),
                         nonZeros[dim1].end(),
                         iv,
    [](const IndVal & iv1, const IndVal & iv2) {
        return iv1.ind < iv2.ind;
    }
                        );
    // itr will point to the end if all indices are less than dim2.
    if (itr == nonZeros[dim1].end())
        return false;
    // itr now points either to dim2 or the next entry after it.
    if ((*itr).ind == dim2) {
        val = &(itr->val);
        return true;
    } else
        return false;
}
} // end namespace SpaMtrix
