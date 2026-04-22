#include <iostream>
// SpaMtrix headers.
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


    // Default finite-difference grid side length is 10 points.
    idx gridLen =5;
    if (nargs > 1){
        gridLen = atoi(args[1]);
    }

    // Create performance timers with millisecond accuracy.
    TickCounter<std::chrono::milliseconds> stopWatch;

    // Create the 5-point Poisson test matrix.
    idx numDoF = gridLen*gridLen;   // Number of degrees of freedom.
    cout << "creating test matrix of size "<<numDoF <<"x" <<numDoF <<"...";

    stopWatch.start();
    MatrixMaker mm(numDoF,numDoF);
    mm.poisson5Point();             // Set the sparsity pattern.
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

    if (gridLen <= 5){ // Print a small grid on screen.
        x.print("x");
    }

    // Print the numerical error magnitude.
    real e = sqrt(errorNorm2(A,x,b));
    cout<< "error norm is " << e << endl;
    cout<< "iterations used " << isol.maxInnerIter << endl;

    // Write the result to a comma-separated text file.
    Writer::writeCSV("out.csv", x , gridLen );
    return 0;
}
