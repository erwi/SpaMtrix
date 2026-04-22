// Copyright (c) 2026 Eero
// SPDX-License-Identifier: MIT
#include <iostream>


// SpaMtrix headers.
#include <spamtrix_matrixmaker.hpp>
#include <spamtrix_ircmatrix.hpp>
#include <spamtrix_vector.hpp>
#include <spamtrix_iterativesolvers.hpp>
#include <spamtrix_diagpreconditioner.hpp>

using std::cout;
using std::endl;
using namespace SpaMtrix;
int main ()
{
  // Example that builds a small sparse matrix and solves it with PCG.
  // A = |3,2|
  //     |2,6|
    
    MatrixMaker mm(2,2);
    mm.addNonZero(0,0,3); 
    mm.addNonZero(0,1,2);
    mm.addNonZero(1,0,2); 
    mm.addNonZero(1,1,6);
    IRCMatrix A = mm.getIRCMatrix();
    cout << "Solving Ax=b, where\nA is:"<<endl;
    A.print();
    
    cout << "b is : "<< endl;
    // Make vectors b and x.
    Vector b(2); b[0] = 2; b[1] = -8;
    Vector x(2);      
    b.print("b");
    
    // Create a diagonal preconditioner.
    DiagPreconditioner M(A);
     // Solve.
    IterativeSolvers isol = IterativeSolvers(10,1e-7);
    bool conv = isol.pcg(A, x, b, M);
    
    std::cout<<"convergence : ";
    if (conv)
      cout << "YES" << endl;
    else
      cout << "NO" << endl;
    
    cout << "maxIter = " << isol.maxIter << endl;
    cout << "toler = " << isol.toler << endl;    
    cout << "solution vector : " << endl;
    x.print("x");
    return 0;
}
