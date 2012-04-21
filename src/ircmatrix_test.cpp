#include <ircmatrix.h>
#include <vector.h>
#include <iostream>
#include <spamtrix_blas.h>
#include <cg.h>
using std::cout;
using std::endl;
int main ()
{
/*!
 * Solves the linear example problem Ax=b from schewchuck's text
 * using the conjugate gradient method.
 */	
    // CREATE 2x2 SPARSE MATRIX A
    idx t[2] = {0,1};
    IRCMatrix A = makeIRCMatrix(t , 2, 1);
    
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
    
    if (conv)
      printf("convergence : YES\n");
    else
      printf("convergence : NO\n");
    
    cout << "maxIter = " << maxIter << endl;
    cout << "toler = " << toler << endl;    
    cout << "solution vector : " << endl;
    x.print();
    return 0;
}