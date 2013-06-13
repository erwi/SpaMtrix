#include <iostream>
#include <cstdlib>
#include <ctime>
#include <spamtrix_matrixmaker.hpp>
#include <spamtrix_ircmatrix.hpp>
#include <spamtrix_vector.hpp>


int main(int nargs, char *args[])
{
    /*!
     * A simple example where a sparse matrix with random sparsity pattern and values is created. \n
     * Tests matrix-vector and matrix-scalar multiplications.
     */
    // DEFAULT MATRIX SIZE IS 5, CAN BE CANGED WITH FIRST COMMAND LINE PARAMETER
    idx testSize = 5;
    if (nargs>1)
        testSize = atoi(args[1]);

    std::cout << "Matrix test size" << testSize << std::endl;

    // CREATE TEST MATRIX WITH RANDOM DATA
    std::cout << "Creating sparse matrix A with random data" << std::endl;
    SpaMtrix::MatrixMaker mm(testSize, testSize);
    srand(time(NULL));
    // FOR EACH ROW AND COLUMN, FILL ~%50 NON_ZEROS
    for (idx r = 0 ; r < testSize ; r++)
        for (idx c = 0 ; c < testSize ; c++){
            if (rand() % 2 ){
                real val = -1.0 + ( (real) (rand() % 1000) ) / 500.0; // RANDOM VALUE IN -1 -> +1 RANGE
                mm.addNonZero(r, c, val);
            }
        }
    SpaMtrix::IRCMatrix Atemp = mm.getIRCMatrix();

    // ASSIGNMENT TEST
    SpaMtrix::IRCMatrix A;
    A = Atemp;
    // DISPLAY CREATED MATRIX AND ITS SPARSITY PATTERN ON SCREEN
    A.print();
    A.spy();

    // MATRIX-VECTOR MULTIPLICATION
    std::cout << "Performing Ax = b, where x is all ones" << std::endl;
    SpaMtrix::Vector x(testSize);
    x = 1.0;
    SpaMtrix::Vector b = A * x;
    b.print("b");

    // MATRIX SCALAR OPERATIONS
    SpaMtrix::IRCMatrix M = -10*A;
    M*=-0.1;
    M.print();



    return 0;
}

