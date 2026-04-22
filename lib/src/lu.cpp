#include <spamtrix_lu.hpp>
#include <spamtrix_vector.hpp>
namespace SpaMtrix {
LU::LU(const IRCMatrix &A):
    numRows(A.getNumRows()) {
  assert(A.getNumRows() == A.getNumCols());
  idx n = A.getNumRows();
#ifdef USES_OPENMP
  #pragma omp sections
    {
        #pragma omp section
        fillFirstColumnL(A, L, n);
        #pragma omp section
        fillFirstRowU(A, U, n);
    }
#else
  fillFirstColumnL(A, L, n);
  fillFirstRowU(A, U, n);
#endif
// Main loop.
  for (idx j = 1; j < n; ++j) {
    // Build L.
    for (idx i = j ; i < n; ++i) {
      real sum(0.0);
      // Accumulate sum = L(i, k) * U(k, j) for k in range 0 to j.
      if (i < L.getNumRows()) {
        for (const auto &nonZero: L.row(i)) {
          const idx k = nonZero.ind;
          if (k >= j) { break; }
          sum += nonZero.val * U.getValue(k, j);
        }
      }
      L.setValue(i, j, A.getValue(i, j) - sum);
    }

    // Build U.
    U.setValue(j, j, 1.0);
    real Ljj = L.getValue(j, j);
    for (idx i = j + 1; i < n; ++i) {
      real sum(0.0);

      for (const auto &nnz : L.row(j)) {
        idx k = nnz.ind;
        if (k >= j) { break; }
        sum += nnz.val * U.getValue(k, i);
      }
      U.setValue(j, i, (A.getValue(j, i) - sum) / Ljj);
    }
  }
}

LU::~LU() { }

void LU::print() {
    std::cout << "L:" << std::endl;
    L.print();
    std::cout << "U:" << std::endl;
    U.print();
}// end void print()

void LU::solve(Vector &x, const Vector &b) const {
    assert(x.getLength() == b.getLength());
    assert(this->numRows == x.getLength());
    Vector y(b.getLength());
    forwardSubstitution(y, b);
    backwardSubstitution(x, y);
}

void LU::forwardSubstitution(Vector &x, const Vector &b) const {
    // Iterate over each row.
    for (idx i = 0 ; i < x.getLength() ; ++i) {
      // Compute the row dot product.
        real sum(0.0);
        for (const auto &nonZero : L.row(i)) {
            sum += nonZero.val * x[nonZero.ind];
        }
        const real diag = L.row(i).back().val; // The diagonal value is the last entry in the row.
        x[i] = (b[i] - sum) / diag;
    }
}

void LU::backwardSubstitution(Vector &x, const Vector &b) const {
    idx n = b.getLength();
    for (idx i = n - 1 ; i < n ; --i) {
        x[i] = b[i];
        real xi = x[i];
        real diag = 0;
        for (const auto &nonZero : U.row(i)) {
            if (nonZero.ind == i) {
                diag = nonZero.val;
            } else {
                xi -= nonZero.val * x[nonZero.ind];
            }
        }
        x[i] = xi / diag;
    }
}
}
