#include <spamtrix_cholesky.hpp>
#include <spamtrix_vector.hpp>

namespace SpaMtrix {
Cholesky::Cholesky(const IRCMatrix &A) {

    for (idx r = 0; r < A.getNumRows(); r++) {
        real Arr = A.sparse_get(r,r);
        // Lower diagonal only, i.e. c < r.
        for (idx c = 0; c <= r; c++) {
            if (r == c) { // Diagonal term.
                real s = 0;
                // Sum of squared values in the row, except for the first row.
                if (r > 0) {
                    for (auto & itr1 : L.row(r)) {
                        s += itr1.val * itr1.val;
                    }
                }

                s = sqrt(Arr - s );

                assert(s > 0.0); // Matrix A is not positive definite.
                L.addNonZero(r, c, s);
            } else { // Off-diagonal.
                real s(0.0);
#ifdef USES_OPENMP
#pragma omp parallel for reduction (+:s)
#endif
                for (idx k = 0; k < c; k++) {
                    real Lrk = L.getValue(r,k);
                    if ( Lrk != 0.0)
                        s+= Lrk*L.getValue(c,k);
                }

                const real a = A.getValue(r, c);
                s =  (a - s);

                if ( s != 0.0 ) { // IF NOT ZERO, ADD TERM TO MATRIX
                    real last= L.row(c).back().val;
                    s /= last;
                    L.addNonZero(r, c, s);
                }
            }
        }
    }
}

Cholesky::~Cholesky() { }

void Cholesky::print() const {
    L.print();
}

void Cholesky::solve(Vector &x, const Vector &b) const {

    Vector y( b.getLength() ); // Temporary vector.
    forwardSubstitution(y,b);
    backwardSubstitution(x,y);
}

void Cholesky::forwardSubstitution(Vector &x, const Vector &b) const {
    // Iterate over each row.
    for (idx i = 0 ; i < x.getLength() ; ++i) {
        // Calculate the sum of the lower-triangular values.
        real sum(0.0);
        auto &cvPairs = L.row(i);

        idx n = (idx) cvPairs.size()-1;
        // Sum all entries below the diagonal.
        for (idx j = 0 ; j < n ; ++j) {
            const idx col = cvPairs[j].ind;
            const real val = cvPairs[j].val;
            sum += x[col] * val;
        }

        // Divide by the diagonal value.
        x[i] = (b[i]-sum) / cvPairs[n].val;
    }
}

void Cholesky::backwardSubstitution(Vector &x, const Vector &b) const {
    idx n = b.getLength();
    // Iterate backward over the rows.
    for (idx i = n-1 ; i < n ; --i) {
        // Perform the sum of L'x, where L' is the transpose of L.
        real xi = b[i];
        for (idx j = i+1 ; j < n ; ++j ) {
            xi -= L.getValue(j,i) * x[j];
        }
        x[i] = xi / L.getValue(i,i);
    }
}
}