#include <ircmatrix.h>
#include <iostream>
using std::cout;
using std::endl;
int main ()
{
	
    idx t[2] = {0,1};
    IRCMatrix A = makeIRCMatrix(t , 2, 1);
    
    // MAKE A = |3,2|
    //		|2,6|
    A.sparse_set(0,0, 3.0);
    A.sparse_set(0,1, 2.0);
    A.sparse_set(1,0, 2.0);
    A.sparse_set(1,1, 6.0);
    
    //A.spy();
    A.print();
    
    
    return 0;
}