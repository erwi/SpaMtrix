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
#include <cholesky.h>
#include <spamtrix_blas.h>
using std::cout;
using std::endl;
int main( int nargs, char *args[] )
{

    // DEFAULT FD DRID SIDE LENGTH IS 10 POINTS
    idx gridLen =10;
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

    cout << " solving Ax = b ..." << endl;
    Cholesky solver(A);
    solver.solve(x,b);
    cout << "OK\n" << endl;

    if (gridLen <= 5) // PRINT SMALL GRID ON SCREEN
        x.print("x");

    // TEST ERROR
    // PRINT NUMERICAL ERROR MAGNITUDE
    real e = sqrt(errorNorm2(A,x,b));
    cout<< "error is " << e << endl;
    
    return 0;

}
