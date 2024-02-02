#include <catch.h>
#include "spamtrix_ircmatrix.hpp"
#include "spamtrix_matrixmaker.hpp"
#include "spamtrix_vector.hpp"
#include "spamtrix_diagpreconditioner.hpp"
#include "spamtrix_iterativesolvers.hpp"
#include "spamtrix_lu.hpp"
#include "spamtrix_cholincpreconditioner.hpp"
#include "spamtrix_luincpreconditioner.hpp"
#include "spamtrix_cholesky.hpp"

void assertEqual(const SpaMtrix::Vector &expected, const SpaMtrix::Vector &actual, double tolerance) {
    REQUIRE(expected.getLength() == actual.getLength());

    for (int i = 0; i < expected.getLength(); i++) {
        double diff = std::abs(expected[i] - actual[i]);
        if (diff > tolerance) {
            std::cerr << diff << " > " << tolerance << " for element " << i << std::endl;
        }
        REQUIRE(expected[i] == Approx(actual[i]).epsilon(tolerance));
    }
}

TEST_CASE("Solve symmetric positive definite 5-point Poisson finite difference system") {
    using namespace SpaMtrix;

    int gridLen = 5;
    idx numDoF = gridLen * gridLen;
    MatrixMaker mm(numDoF,numDoF);
    mm.poisson5Point();

    IRCMatrix A = mm.getIRCMatrix();

    Vector b(numDoF);
    b.setAllValuesTo(1.);

    // expected solution values obtained using Octave A \ b operation
    Vector xpected(25);
    xpected[0] = 0.951923076923077;
    xpected[1] = 1.403846153846153;
    xpected[2] = 1.538461538461538;
    xpected[3] = 1.403846153846154;
    xpected[4] = 0.951923076923077;
    xpected[5] = 1.403846153846153;
    xpected[6] = 2.125000000000000;
    xpected[7] = 2.346153846153846;
    xpected[8] = 2.125000000000000;
    xpected[9] = 1.403846153846154;
    xpected[10] = 1.538461538461538;
    xpected[11] = 2.346153846153845;
    xpected[12] = 2.596153846153845;
    xpected[13] = 2.346153846153845;
    xpected[14] = 1.538461538461538;
    xpected[15] = 1.403846153846153;
    xpected[16] = 2.124999999999999;
    xpected[17] = 2.346153846153845;
    xpected[18] = 2.124999999999999;
    xpected[19] = 1.403846153846154;
    xpected[20] = 0.951923076923077;
    xpected[21] = 1.403846153846153;
    xpected[22] = 1.538461538461538;
    xpected[23] = 1.403846153846153;
    xpected[24] = 0.951923076923077;

    const double eps = 1e-15;

    SECTION("Solve using LU decomposition") {
        Vector x(b.getLength());
        LU lu(A);

        lu.solve(x, b);

        assertEqual(xpected, x, 10 * eps);
    }

    SECTION("Solve using Cholesky decomposition") {
        Vector x(b.getLength());
        Cholesky cholesky(A);

        cholesky.solve(x, b);

        assertEqual(xpected, x, 10 * eps);
    }

    SECTION("Solve using PCG and diagonal preconditioner") {
        IterativeSolvers solver(numDoF, eps);
        DiagPreconditioner preconditioner(A);

        Vector x(b.getLength());
        bool converged = solver.pcg(A, x, b, preconditioner);

        REQUIRE(converged);
        REQUIRE(solver.toler <= eps);
        assertEqual(xpected, x, 10 * eps);
    }

    SECTION("Solve using PCG and incomplete Cholesky preconditioner") {
        IterativeSolvers solver(numDoF, eps);

        CholIncPreconditioner preconditioner(A);

        Vector x(b.getLength());
        bool converged = solver.pcg(A, x, b, preconditioner);

        REQUIRE(converged);
        REQUIRE(solver.toler <= eps);
        assertEqual(xpected, x, 10 * eps);
    }

    SECTION("Solve using PCG and incomplete LU preconditioner") {
        IterativeSolvers solver(numDoF, eps);
        LUIncPreconditioner preconditioner(A);
        Vector x(b.getLength());
        bool converged = solver.pcg(A, x, b, preconditioner);

        REQUIRE(converged);
        REQUIRE(solver.toler <= eps);
        assertEqual(xpected, x, 10 * eps);
    }

    SECTION("Solve using GMRES and diagonal preconditioner") {
        IterativeSolvers solver(numDoF, 100, eps);
        DiagPreconditioner preconditioner(A);

        Vector x(b.getLength());
        bool converged = solver.gmres(A, x, b, preconditioner);

        REQUIRE(converged);
        REQUIRE(solver.toler <= eps);
        assertEqual(xpected, x, 10 * eps);
    }

    SECTION("Solve using GMRES and incomplete Cholesky preconditioner") {
        IterativeSolvers solver(numDoF, 100, eps);
        CholIncPreconditioner preconditioner(A);

        Vector x(b.getLength());
        bool converged = solver.gmres(A, x, b, preconditioner);

        REQUIRE(converged);
        REQUIRE(solver.toler <= eps);
        assertEqual(xpected, x, 10 * eps);
    }

    SECTION("Solve using GMRES and incomplete LU preconditioner") {
        IterativeSolvers solver(numDoF, 100, eps);
        LUIncPreconditioner preconditioner(A);

        Vector x(b.getLength());
        bool converged = solver.gmres(A, x, b, preconditioner);

        REQUIRE(converged);
        REQUIRE(solver.toler <= eps);
        assertEqual(xpected, x, 10 * eps);
    }
}
