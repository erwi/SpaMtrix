#include <spamtrix_cholincpreconditioner.hpp>
namespace SpaMtrix {
CholIncPreconditioner::CholIncPreconditioner(const IRCMatrix &A) {
    /*!
    * Incoplete cholesky matrix constructor.
    */
    L.numDim1 = A.getNumRows();
    L.numDim2 = A.getNumCols();
    // FOR EACH ROW
    for (idx r = 0 ; r < A.getNumRows() ; r++) {
        L.nonZeros.push_back(std::vector<IndVal>());  // ADD NEW EMPTY ROW
        // FOR EACH COLUMN, LOWERD DIAGONAL ONLY, i.e. c < r
        for (idx c = 0 ; c <= r ; c++) {
            if (r == c) { // DIAGONAL
                real s = 0;
                // SUM OVER ALL VALUES SQUARED TO LEFT ON THIS ROW
                std::vector< IndVal>::iterator itr1 = L.nonZeros[r].begin();
                std::vector< IndVal>::iterator itr2 = L.nonZeros[r].end();
                for (; itr1 != itr2 ; itr1++)
                    s += itr1->val * itr1->val;
                s = sqrt(A.sparse_get(r, r) - s);
#ifdef DEBUG
                if (s <= 0.0) {
                    std::cerr << "error in " << __func__ << "matrix A is not positive definite - bye!" << std::endl;
                    exit(1);
                }
#endif
                L.nonZeros[r].push_back(IndVal(c, s));// std::move could be used here?
            } else { // OFF-DIAGONAL
                real a(0.0);
                if (!A.isNonZero(r, c, a))
                    continue;
                real s = 0;
                for (idx k = 0 ; k < c; k++)
                    s += L.getValue(r, k) * L.getValue(c, k); // OPTIMISE ITERATION OVER NONZEROS IN ROW r
                s = (a - s);
                if (s != 0.0) { // IF NOT ZERO, ADD TERM TO MATRIX
                    real last = L.nonZeros[c].back().val;
                    s /= last;
                    L.nonZeros[r].push_back(IndVal(c, s));
                }
            }
        }//end for columns
    }// end for rows
}
CholIncPreconditioner::~CholIncPreconditioner() { }

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
        std::vector<IndVal>::const_iterator itr1 = L.nonZeros[i].begin();
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
