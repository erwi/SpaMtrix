#include <iostream>
#include <string>

#include <reader.h>
#include <writer.h>
#include <matrixmaker.h>
#include <ircmatrix.h>
/*! Example that reads/writes matrices from/to files.
    1. First create a matrix. 
    2. Write it into a file.
    3. Read the file into another matrix.
 */
int main(int nargs, char* args[] )
{
    // 1. CREATE TEST MATRIX
    SpaMtrix::MatrixMaker mm(25,25);
    mm.poisson5Point();
    SpaMtrix::IRCMatrix A = mm.getIRCMatrix();
    
    A.spy();
    
    // 2. WRITE TEST MATRIX ONTO FILE
    std::string filename("testoutput.txt");
    SpaMtrix::Writer::writeMatrixMarket(filename,A);
    A.clear();
    
    // 3. READ TEST MATRIX FROM FILE
    
    
    
    return 0;
}