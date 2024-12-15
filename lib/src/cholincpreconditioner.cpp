#include <spamtrix_cholincpreconditioner.hpp>
#include <spamtrix_exception.hpp>

namespace SpaMtrix {
CholIncPreconditioner::CholIncPreconditioner(const IRCMatrix &A) : L() {
  for (idx r = 0 ; r < A.getNumRows() ; r++) {
    // for each column, lower diagonal only, i.e. c < r
    for (idx c = 0 ; c <= r ; c++) {
      if (r == c) { // DIAGONAL
        real s = 0;
        // Sum of squared row values
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
    /*!
         x[i] =  ( b[i] - sum( A[i,j]x[j] ) ) / A[i,i]
     */
    // FOR EACH ROW, STARTING FROM FIRST
    for (idx i = 0 ; i < x.getLength() ; ++i) {
        // FORM SUM OVER ALL LOWER DIAGONAL MATRIX VALUES
        real sum(0.0);
        auto itr1 = L.row(i).begin();
        idx col = itr1->ind;
        // WHILE LOWER DIAGONAL NONZEROS ONLY
        // UPPER DIAGONAL VALUES TAKEN CARE OF IN back-substitution
        while (col < i) {
            sum += x[col] * (itr1->val);
            itr1++;
            col = itr1->ind;
        }
        x[i] = (b[i] - sum) / itr1->val;
    }
}

void CholIncPreconditioner::backwardSubstitution(Vector &x, const Vector &b) const {
    /*!
      Performs back substitution of transposed lower triangular matrix.

      x[i] = ( b[i] - sum(L'x[i+1:end]) ) / L(i,i)
    */
    idx n = b.getLength();
    // FOR EACH ROW, STARTING FROM LAST, COUNTING BACKWARDS
    // NOTE: i < n FOR UNSIGNED INTS IS EQUIVALENT TO
    // i >= 0 FOR SIGNED INTS
    for (idx i = n - 1 ; i < n ; --i) {
        x[i] = b[i];
        // PERFORM SUM OF L'x, WHERE L' IS TRANSPOSE OF L
        // THIS USES MATRIX SEARCH AND SHOULD BE OPTIMISED
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
