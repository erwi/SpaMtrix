#include <spamtrix_lu.hpp>
#include <spamtrix_vector.hpp>
namespace SpaMtrix {
LU::LU(const IRCMatrix &A):
    numRows(A.getNumRows()),
    L(numRows, numRows),
    U(numRows, numRows) {
#ifdef DEBUG
    assert(A.getNumRows() == A.getNumCols());
#endif
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
//==============================
//          MAIN LOOP
//==============================
    for (idx j = 1; j < n ; ++j) {
        // FOR ROWS i
        for (idx i = j ; i < n ; ++i) {
            real sum(0.0);
#ifdef USES_OPENMP
            #pragma omp parallel for reduction(+:sum)
#endif
            for (idx k = 0; k < j ; ++k) {
                real Lik = L.getValue(i, k);
                if (Lik != 0.0) {
                    sum += Lik * U.getValue(k, j);
                }
            }//end for k
            L.setValue(i, j, A.getValue(i, j) - sum);
        }// end for i
        // MAKE U
        U.setValue(j, j, 1.0);
        real Ljj = L.getValue(j, j);
        for (idx i = j + 1; i < n; ++i) {
            real sum(0.0);
#ifdef USES_OPENMP
            #pragma omp parallel for reduction(+:sum)
#endif
            for (idx k = 0 ; k < j ; ++k) {
                real Ljk = L.getValue(j, k);
                if (Ljk != 0.0) {
                    sum += Ljk * U.getValue(k, i);
                }
            }//end for k
            U.setValue(j, i, (A.getValue(j, i) - sum) / Ljj);
        }// end for i
    }// end for j
}// end constructor

LU::~LU() { }

void LU::print() {
    /*!
    * Debug printout to stdout of L and U Matrices.
    */
    std::cout << "L:" << std::endl;
    L.print();
    std::cout << "U:" << std::endl;
    U.print();
}// end void print()

void LU::solve(Vector &x, const Vector &b) const {
    /*!
    * Solves Ax=b using forward/backward substitution. \n
    * Ax  = b \n
    * L(Uy) = b \n
    * Lx = y \n
    */
    Vector y(b.getLength());
    forwardSubstitution(y, b);
    backwardSubstitution(x, y);
}// end void solve

void LU::forwardSubstitution(Vector &x, const Vector &b) const {
    // FOR EACH ROW
    for (idx i = 0 ; i < x.getLength() ; ++i) {
        // FORM SUM
        real sum(0.0);
        // REFERENCE TO ALL NON-ZEROS IN ROW i
        const std::vector<IndVal> &cvPairs = L.nonZeros[i];
        idx n = (idx) cvPairs.size() - 1; // ROW LENGTH, NOT INCLUDING DIAGONAL
#ifdef USES_OPENMP
        #pragma omp parallel for reduction(+:sum)
#endif
        for (idx j = 0 ; j < n ; ++j) {
            const idx col = cvPairs[j].ind;
            const real val = cvPairs[j].val;
            sum += x[col] * val;
        }
        x[i] = (b[i] - sum) / cvPairs[n].val; // DIVIDE BY DIAGONAL VALUE
    }
}// end void forwardSubstitution

void LU::backwardSubstitution(Vector &x, const Vector &b) const {
    idx n = b.getLength();
    for (idx i = n - 1 ; i < n ; --i) {
        x[i] = b[i];
        real xi = x[i];
        const std::vector<IndVal> &cvPairs = U.nonZeros[i];
        const idx n = cvPairs.size();
        for (idx j = 1 ; j < n ; ++j) {
            const idx col = cvPairs[j].ind;
            const real val = cvPairs[j].val;
            xi -= val * x[col];
        }
        x[i] = xi / cvPairs[0].val;
    }// end for i
}// end void backwardSubstitution
} // end namespace SpaMtrix
