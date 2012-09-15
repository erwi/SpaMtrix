#include <iostream>


// SpaMtrix HEADERS
#include <matrixmaker.h>
#include <ircmatrix.h>
#include <vector.h>
#include <cg.h>
#include <diagpreconditioner.h>

using std::cout;
using std::endl;
int main ()
{
/*!
 * Example/test showing how to create a sparse matrix using the MatrixMaker.
 * 
 * A small 2x2 matrix problem Ax=b is constructed and solver using the 
 * preconditioned conjugate gradient method.
 * 
 * A = [[3,2];[2,6]] , b = [2,-8] is used.
 * 
 */	
    // CREATE SPARSE MATRIX    
    //      A = |3,2|
    //          |2,6|
    
    MatrixMaker mm(2,2);
    mm.addNonZero(0,0,3); 
    mm.addNonZero(0,1,2);
    mm.addNonZero(1,0,2); 
    mm.addNonZero(1,1,6);
    IRCMatrix A = mm.getIRCMatrix();
    cout << "Solving Ax=b, where\nA is:"<<endl;
    A.print();
    
    cout << "b is : "<< endl;
    // MAKE VECTORS b AND x
    Vector b(2); b[0] = 2; b[1] = -8;
    Vector x(2);      
    b.print("b");
    
    // SOLVE
    idx maxIter = 10; // MAXIMUM NUMBER OF ITERATION USED
    real toler(1e-7); // REQUIRED NUMERICAL ACCURACY 
    
    // CREATE DIAGONAL PRECONDITIONER MATRIX
    DiagPreconditioner M(A);
    
    bool conv = cg(A, x, b, M, maxIter, toler );
    
    std::cout<<"convergence : ";
    if (conv)
      cout << "YES" << endl;
    else
      cout << "NO" << endl;
    
    cout << "maxIter = " << maxIter << endl;
    cout << "toler = " << toler << endl;    
    cout << "solution vector : " << endl;
    x.print("x");
    return 0;
}