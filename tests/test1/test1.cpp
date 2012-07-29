#include <ircmatrix.h>
#include <vector.h>
#include <iostream>
#include <spamtrix_blas.h>
#include <cg.h>
#include <matrixmaker.h>
using std::cout;
using std::endl;
int main ()
{
/*!
 * Solves the linear example problem Ax=b from schewchuck's text
 * using the conjugate gradient method.
 */	
    // CREATE EMPTY SPARSE MATRIX    
    MatrixMaker mm;
    mm.setMatrixSize(2,2);
    mm.addNonZero(0,0); mm.addNonZero(0,1);
    mm.addNonZero(1,0); mm.addNonZero(1,1);
    IRCMatrix A = mm.getIRCMatrix();
    mm.clear();
    
    // SET  A = |3,2|
    //		|2,6|
    A.sparse_set(0,0, 3.0);
    A.sparse_set(0,1, 2.0);
    A.sparse_set(1,0, 2.0);
    A.sparse_set(1,1, 6.0);
    A.print();
    
    // MAKE VECTORS b AND x
    Vector b(2); b[0] = 2; b[1] = -8;
    Vector x(2);      
    
    // SOLVE
    idx maxIter = 10;
    real toler(1e-7);
    bool conv = cg(A, b, x, maxIter, toler );
    
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