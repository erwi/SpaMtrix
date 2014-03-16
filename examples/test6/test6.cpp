/*!
  Solution of system Ax = b.

  Creates a 5 point poisson (2D finite differences) test matrix A.

  The FD grid is assumed square, with equal side lengths n, and
  grid spacing of unity.

  The boundary conditions are all zero Dirchlet nodes and not solved for.

  The right hand side vector b is set to all ones.

  The result vector x is writen out in a comma separated values text file 'out.txt'
  that can easily be loaded and visualised with a spreadsheet/MATLAB/OCTAVE etc. program.
*/

#include <iostream>
// SpaMtrix HEADERS
#include <spamtrix_matrixmaker.hpp>
#include <spamtrix_ircmatrix.hpp>
#include <spamtrix_vector.hpp>
#include <spamtrix_blas.hpp>
#include <spamtrix_iterativesolvers.hpp>
#include <spamtrix_writer.hpp>
#include <spamtrix_cholincpreconditioner.hpp>
#include <spamtrix_diagpreconditioner.hpp>
#include <spamtrix_tickcounter.hpp>
using std::cout;
using std::endl;
using namespace SpaMtrix;
int main( int nargs, char *args[] )
{
#ifdef USES_OPENMP
omp_set_num_threads(0);
#endif


    // DEFAULT FD GRID SIDE LENGTH IS 10 POINTS
    idx gridLen =5;
    if (nargs > 1){
        gridLen = atoi(args[1]);
    }

    // CREATE PERFORMANCE TIMES WITH MILLISECOND ACCURACY
    TickCounter<std::chrono::milliseconds> stopWatch;

    // CREATE 5-POINT-POISSON TEST MATRIX
    idx numDoF = gridLen*gridLen;   // NUMBER OF DEGREES OF FREEDOM OF SYSTEM
    cout << "creating test matrix of size "<<numDoF <<"x" <<numDoF <<"...";

    stopWatch.start();
    MatrixMaker mm(numDoF,numDoF);
    mm.poisson5Point();             // SETS SPARSITY PATTERN
    //IRCMatrix A = mm.getIRCMatrix();
    IRCMatrix A;
    mm.makeSparseMatrix(A);
    cout << "OK, elapsed: " << stopWatch.getElapsed() << endl;

    // CREATE VECTORS FOR SYSTEM OF EQUATIONS Ax = b
    Vector x(numDoF);
    Vector b(numDoF);
    b = 1.0;


    stopWatch.reset();
    cout << "making preconditioner ...";
    DiagPreconditioner M(A);
    cout << "OK, time elapsed " << stopWatch.getElapsed() << "ms" << endl;
    cout << " solving Ax = b ..." << endl;

    stopWatch.reset();
    IterativeSolvers isol(numDoF, gridLen, 1e-7);
    isol.gmres(A, x, b, M);
    cout << "OK, solved in " << stopWatch.getElapsed() << "ms" << endl;

    if (gridLen <= 5){ // PRINT SMALL GRID ON SCREEN
        x.print("x");
    }

    // TEST ERROR
    // PRINT NUMERICAL ERROR MAGNITUDE
    real e = sqrt(errorNorm2(A,x,b));
    cout<< "error norm is " << e << endl;
    cout<< "iterations used " << isol.maxInnerIter << endl;

    // WRITE RESULT IN A COMMA SEPARATED TEXT FILE
    // ROWS AND COLUMNS ARE ORDERED ACCORDING TO THE
    // FD GRID USED FOR THE CALCULATION
    Writer::writeCSV("out.csv", x , gridLen );
    return 0;
}
