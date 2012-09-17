/*!
  Creating a 5 point poisson (2D finite differences) test matrix.
  The FD grid is assumed square, with equal side lengths n, and
  grid spacing of unity.
*/

#include <iostream>
#include <omp.h>
// SpaMtrix headers
#include <matrixmaker.h>
#include <ircmatrix.h>
#include <vector.h>
#include <spamtrix_blas.h>
#include <gmres.h>
#include <writer.h>
#include <cholincpreconditioner.h>
#include <diagpreconditioner.h>
#include <sorpreconditioner.h>
#include <densematrix.h>
#include <tickcounter.h>	// PERFORMANCE TIMER
using std::cout;
using std::endl;
int main( int nargs, char *args[] )
{
omp_set_num_threads(0);
    // DEFAULT FD GRID SIDE LENGTH IS 10 POINTS
    idx gridLen =5;
    if (nargs > 1)
    {
        gridLen = atoi(args[1]);
    }
 
    // CREATE 5-POINT-POISSON TEST MATRIX
    idx numDoF = gridLen*gridLen;   // NUMBER OF DEGREES OF FREEDOM OF SYSTEM
    cout << "creating test matrix of size "<<numDoF <<"x" <<numDoF <<"...";
    MatrixMaker mm(numDoF,numDoF);
    mm.poisson5Point();             // SETS SPARSITY PATTERN
    IRCMatrix A = mm.getIRCMatrix();
    cout << "OK" << endl;
    
    // CREATE VECTORS FOR SYSTEM OF EQUATIONS Ax = b
    Vector x(numDoF);
    Vector b(numDoF);
    b[0] = 1.0;

    // CREATE PERFORMANCE TIMES WITH MILLISECOND ACCURACY
    TickCounter<std::chrono::milliseconds> stopWatch;
    
    stopWatch.start();
    cout << "making preconditioner ...";
    //CholIncPreconditioner M(A);
    DiagPreconditioner M(A);
    //SORPreconditioner M(A,20);
    cout << "OK, time elapsed " << stopWatch.getElapsed() << "ms" << endl;
    
    cout << " solving Ax = b ..." << endl;

    idx maxIter(numDoF);
    idx innerIter(gridLen);
    real toler(1e-6);
    DenseMatrix H(innerIter+1, innerIter);
    
    stopWatch.reset();
    gmres(A, x, b, M, H, innerIter, maxIter, toler);
    
    cout << "OK, solved in " << stopWatch.getElapsed() << "ms" << endl;

    if (gridLen <= 5) // PRINT SMALL GRID ON SCREEN
        x.print("x");

    // TEST ERROR
    // PRINT NUMERICAL ERROR MAGNITUDE
    real e = sqrt(errorNorm2(A,x,b));
    cout<< "error norm is " << e << endl;
    cout<< "iterations used " << maxIter << endl;
    
    // WRITE RESULT IN A COMMA SEPARATED TEXT FILE
    // ROWS AND COLUMNS ARE ORDERED ACCORDING TO THE
    // FD GRID USED FOR THE CALCULATION
    Writer w;
    w.writeCSV("out.csv", x , gridLen );
    
    return 0;

}
