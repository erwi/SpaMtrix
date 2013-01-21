/*!
  Creating a 5 point poisson (2D finite differences) test matrix.
  The FD grid is assumed square, with equal side lengths n, and
  grid spacing of unity.
*/

#include <iostream>

// SpaMtrix headers
#include <spamtrix_matrixmaker.hpp>
#include <spamtrix_ircmatrix.hpp>
#include <spamtrix_vector.hpp>
#include <spamtrix_cholesky.hpp>
#include <spamtrix_blas.hpp>
#include <spamtrix_writer.hpp>
using std::cout;
using std::endl;
using namespace SpaMtrix;
int main( int nargs, char *args[] )
{

    // DEFAULT FD DRID SIDE LENGTH IS 10 POINTS
    idx gridLen =10;
    if (nargs > 1)
    {
        gridLen = atoi(args[1]);
    }
    cout << "Creating 5-point poisson test matrix A..." << endl;
    // CREATE 5-POINT-POISSON TEST MATRIX
    idx numDoF = gridLen*gridLen;   // NUMBER OF DEGREES OF FREEDOM OF SYSTEM
    MatrixMaker mm(numDoF,numDoF);
    mm.poisson5Point();             // SETS SPARSITY PATTERN
    IRCMatrix A = mm.getIRCMatrix();
    cout << "Matrix size is : " << numDoF << "x" << numDoF << endl;
    // CREATE VECTORS FOR SYSTEM OF EQUATIONS Ax = b
    Vector x(numDoF);
    Vector b(numDoF);
    b[0] = 1.0;

    cout << "solving Ax = b using Cholesky decomposition...";
    Cholesky solver(A);
    solver.solve(x,b);
    cout << "OK" << endl;

    if (gridLen <= 5) // PRINT SMALL GRID ON SCREEN
        x.print("x");

    // TEST ERROR
    // PRINT NUMERICAL ERROR MAGNITUDE
    real e = sqrt(errorNorm2(A,x,b));
    cout<< "error is " << e << endl;
    
    // WRITE RESULT IN A COMMA SEPARATED TEXT FILE
    // ROWS AND COLUMNS ARE ORDERED ACCORDING TO THE
    // FD GRID USED FOR THE CALCULATION
    Writer w;
    w.writeCSV("out.csv", x , gridLen );
    
    return 0;

}
