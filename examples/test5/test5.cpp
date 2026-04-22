// Copyright (c) 2026 Eero
// SPDX-License-Identifier: MIT
#include <iostream>

// SpaMtrix headers.
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

    // Default finite-difference grid side length is 10 points.
    idx gridLen =10;
    if (nargs > 1)
    {
        gridLen = atoi(args[1]);
    }
    cout << "Creating 5-point poisson test matrix A..." << endl;
    // Create the 5-point Poisson test matrix.
    idx numDoF = gridLen*gridLen;   // Number of degrees of freedom.
    MatrixMaker mm(numDoF,numDoF);
    mm.poisson5Point();             // Set the sparsity pattern.
    IRCMatrix A = mm.getIRCMatrix();
    cout << "Matrix size is : " << numDoF << "x" << numDoF << endl;
    // Create vectors for the system of equations Ax = b.
    Vector x(numDoF);
    Vector b(numDoF);
    b[0] = 1.0;

    cout << "solving Ax = b using Cholesky decomposition...";
    Cholesky solver(A);
    solver.solve(x,b);
    cout << "OK" << endl;

    if (gridLen <= 5) // Print a small grid on screen.
        x.print("x");

    // Print the numerical error magnitude.
    real e = sqrt(errorNorm2(A,x,b));
    cout<< "error is " << e << endl;
    
    // Write the result to a comma-separated text file.
    Writer::writeCSV("out.csv", x , gridLen );
    
    return 0;

}
