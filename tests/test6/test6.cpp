/*!
  Creating a 5 point poisson (2D finite differences) test matrix.
  The FD grid is assumed square, with equal side lengths n, and
  grid spacing of unity.
*/

#include <iostream>

// SpaMtrix headers
#include <matrixmaker.h>
#include <ircmatrix.h>
#include <vector.h>
#include <spamtrix_blas.h>
#include <gmres.h>
#include <writer.h>
#include <cholincpreconditioner.h>
#include <densematrix.h>
using std::cout;
using std::endl;
int main( int nargs, char *args[] )
{

    // DEFAULT FD GRID SIDE LENGTH IS 10 POINTS
    idx gridLen =5;
    if (nargs > 1)
    {
        gridLen = atoi(args[1]);
    }
 
    // CREATE 5-POINT-POISSON TEST MATRIX
    idx numDoF = gridLen*gridLen;   // NUMBER OF DEGREES OF FREEDOM OF SYSTEM
    MatrixMaker mm(numDoF,numDoF);
    mm.poisson5Point();             // SETS SPARSITY PATTERN
    IRCMatrix A = mm.getIRCMatrix();

    // CREATE VECTORS FOR SYSTEM OF EQUATIONS Ax = b
    Vector x(numDoF);
    Vector b(numDoF);
    b[0] = 1.0;

    CholIncPreconditioner M(A);

    cout << " solving Ax = b ..." << endl;

    idx maxIter(20);
    idx innerIter(10);
    real toler(1e-6);
    DenseMatrix H(numDoF+1, numDoF);
    bool converged =
    gmres(A, x, b, M, H, innerIter, maxIter, toler);

    cout << "OK\n" << endl;

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
