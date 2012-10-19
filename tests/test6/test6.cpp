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
#include <omp.h>
// SpaMtrix HEADERS
#include <matrixmaker.h>
#include <ircmatrix.h>
#include <vector.h>
#include <spamtrix_blas.h>
#include <iterativesolvers.h>
#include <writer.h>
#include <cholincpreconditioner.h>
#include <diagpreconditioner.h>
#include <sorpreconditioner.h>
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
 
    // CREATE PERFORMANCE TIMES WITH MILLISECOND ACCURACY
    TickCounter<std::chrono::milliseconds> stopWatch;
 
 
    // CREATE 5-POINT-POISSON TEST MATRIX
    idx numDoF = gridLen*gridLen;   // NUMBER OF DEGREES OF FREEDOM OF SYSTEM
    cout << "creating test matrix of size "<<numDoF <<"x" <<numDoF <<"...";
    
    stopWatch.start();
    MatrixMaker mm(numDoF,numDoF);
    mm.poisson5Point();             // SETS SPARSITY PATTERN
    IRCMatrix A = mm.getIRCMatrix();
    cout << "OK, elapsed: " << stopWatch.getElapsed() << endl;
    
    // CREATE VECTORS FOR SYSTEM OF EQUATIONS Ax = b
    Vector x(numDoF);
    Vector b(numDoF);
    b = 1.0;

    
    stopWatch.reset();
    cout << "making preconditioner ...";
    //CholIncPreconditioner M(A);
    DiagPreconditioner M(A);
    
    cout << "OK, time elapsed " << stopWatch.getElapsed() << "ms" << endl;
    
    cout << " solving Ax = b ..." << endl;

    idx maxIter(numDoF);
    idx innerIter(gridLen);
    real toler(1e-6);
    stopWatch.reset();
    IterativeSolvers::gmres(A, x, b, M, maxIter, innerIter, toler);
    
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
