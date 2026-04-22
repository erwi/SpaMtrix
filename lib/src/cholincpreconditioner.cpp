#include <spamtrix_cholincpreconditioner.hpp>
#include <spamtrix_exception.hpp>

namespace SpaMtrix {
CholIncPreconditioner::CholIncPreconditioner(const IRCMatrix &A) : L() {
  for (idx r = 0 ; r < A.getNumRows() ; r++) {
    // for each column, lower diagonal only, i.e. c < r
    for (idx c = 0 ; c <= r ; c++) {
      if (r == c) { // DIAGONAL
        real s = 0;
        // Sum of squared row values.
        if (r != 0) {
          for (auto &nonZero: L.row(r)) {
            s += nonZero.val * nonZero.val;
          }
        }

        s = sqrt(A.sparse_get(r, r) - s);
        if (s <= 0.0) {
          throw SpaMtrixException("error in " + std::string(ERROR_LOCATION) + "matrix A is not positive definite.");
        }

        L.addNonZero(r, c, s);
      } else { // OFF-DIAGONAL
        real a(0.0);
        if (!A.isNonZero(r, c, a)) {
          continue;
        }
        real s = (a - s); // TODO: Check this, s appears to be uninitialised here!!
        if (s != 0.0) { // IF NOT ZERO, ADD TERM TO MATRIX
          real last = L.row(c).back().val;
          s /= last;
          L.addNonZero(r, c, s);
        }
      }
    }
  }
}

CholIncPreconditioner::~CholIncPreconditioner() = default;

void CholIncPreconditioner::solveMxb(Vector &x, const Vector &b) const {
    Vector y(b.getLength());  // TEMPORARY VECTOR
    forwardSubstitution(y, b);
    backwardSubstitution(x, y);
}
void CholIncPreconditioner::forwardSubstitution(Vector &x, const Vector &b) const {
    // Iterate over each row from the beginning.
    for (idx i = 0 ; i < x.getLength() ; ++i) {
      // Form the sum over all lower-diagonal matrix values.
        real sum(0.0);
        auto itr1 = L.row(i).begin();
        idx col = itr1->ind;
        // Stop at the diagonal; upper-diagonal values are handled in back substitution.
        while (col < i) {
            sum += x[col] * (itr1->val);
            itr1++;
            col = itr1->ind;
        }
        x[i] = (b[i] - sum) / itr1->val;
    }
}

void CholIncPreconditioner::backwardSubstitution(Vector &x, const Vector &b) const {
    idx n = b.getLength();
    // Iterate backward over the rows.
    for (idx i = n - 1 ; i < n ; --i) {
        x[i] = b[i];
      // Perform the sum of L'x, where L' is the transpose of L.
        for (idx j = n - 1 ; j > i ; --j) {
            x[i] -= L.getValue(j, i) * x[j];
        }
        x[i] /= L.getValue(i, i);
    }// end for each row
}

void CholIncPreconditioner::print() const {
    L.print();
}
} // end namespace SpaMtrix
