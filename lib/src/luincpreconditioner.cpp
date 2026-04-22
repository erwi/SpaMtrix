// Copyright (c) 2026 Eero
// SPDX-License-Identifier: MIT
#include <spamtrix_luincpreconditioner.hpp>
namespace SpaMtrix {
LUIncPreconditioner::~LUIncPreconditioner() {}

LUIncPreconditioner::LUIncPreconditioner(const IRCMatrix &A) : M(A) {
    const idx n = A.getNumRows();

    for (idx i = 1 ; i < n ; ++i) {
        auto k_itr = M.row(i).begin();
        auto rowEnd = M.row(i).end();
        // Iterate over the lower-triangular entries in row i.
        while (k_itr->ind < i ) { // a[i, k]
            const idx k = k_itr->ind;
            k_itr->val /= M.getValue(k,k);
            // Iterate over the entries to the right of k.
            auto j_itr = k_itr + 1;
            while ( j_itr != rowEnd ) { // a[i, j]
                const idx j = j_itr->ind;
                j_itr->val -= k_itr->val * M.getValue(k,j); // TODO: Optimise row access for non-zeros.
                j_itr++;
            }
            k_itr ++;
        }
    }
}

void LUIncPreconditioner::solveMxb(Vector &x, const Vector &b) const {
    assert(x.getLength() == b.getLength() );
    Vector y( b.getLength());
    forwardSubstitution(y, b);
    backwardSubstitution(x,y);
}

void LUIncPreconditioner::forwardSubstitution(Vector &x, const Vector &b) const {
    const idx n = x.getLength();
    // Iterate over each row.
    for(idx i = 0 ; i < n ; ++i) {
        // Calculate the sum of all lower-diagonal values.
        real sum(0.0);
        auto col_itr = M.row(i).begin();
        while (col_itr->ind < i) {
            sum += x[col_itr->ind] * (col_itr->val);
            ++col_itr;
        }
        x[i] = b[i]-sum;
    }
}

void LUIncPreconditioner::backwardSubstitution(Vector &x, const Vector &b) const {
    const idx n = b.getLength();

    for (idx i = n-1; i < n; --i) {
        x[i] = b[i];
        auto col_itr = M.row(i).rbegin();

        while ( col_itr->ind > i) {
            x[i] -= col_itr->val * x[col_itr->ind];
            ++col_itr;
        }
        x[i] /= col_itr->val; // finally, division by diagonal
    }
}

void LUIncPreconditioner::print() const {
    M.print();
}
} // end namespace
