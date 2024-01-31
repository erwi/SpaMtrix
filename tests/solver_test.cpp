#include <catch.h>
#include "spamtrix_ircmatrix.hpp"
#include "spamtrix_matrixmaker.hpp"
#include "spamtrix_vector.hpp"
#include "spamtrix_diagpreconditioner.hpp"
#include "spamtrix_iterativesolvers.hpp"

TEST_CASE("PCG solver test") {
    // Solve Ax=b using PCG when
    // A = |3, 2|
    //     |2, 6|
    //
    // b = [2, -8]
    //
    // and expect x = [2, -2]

    using namespace SpaMtrix;

    MatrixMaker mm(2,2);
    mm.addNonZero(0,0,3);
    mm.addNonZero(0,1,2);
    mm.addNonZero(1,0,2);
    mm.addNonZero(1,1,6);
    IRCMatrix A = mm.getIRCMatrix();

    Vector b(2); b[0] = 2; b[1] = -8;
    Vector x(2);

    // Solve using pcg and using Diagonal preconditioner
    DiagPreconditioner M(A);
    IterativeSolvers isol = IterativeSolvers(10,1e-7);
    bool converged = isol.pcg(A, x, b, M);

    REQUIRE(converged);
    REQUIRE(x[0] == Approx(2).epsilon(1e-7));
    REQUIRE(x[1] == Approx(-2).epsilon(1e-7));
}
