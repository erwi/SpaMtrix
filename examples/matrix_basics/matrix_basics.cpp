#include <iostream>
#include <cstdlib>
#include <ctime>
#include <spamtrix_matrixmaker.hpp>
#include <spamtrix_ircmatrix.hpp>
#include <spamtrix_vector.hpp>


int main(int nargs, char *args[])
{
    // A simple example that creates a sparse matrix with random values.
    // It also tests matrix-vector and matrix-scalar multiplications.
    // Default matrix size is 5 and can be changed with the first command-line parameter.
    idx testSize = 5;
    if (nargs>1)
        testSize = atoi(args[1]);

    std::cout << "Matrix test size" << testSize << std::endl;

    // Create a test matrix with random data.
    std::cout << "Creating sparse matrix A with random data" << std::endl;
    SpaMtrix::MatrixMaker mm(testSize, testSize);
    srand(time(NULL));
    // Fill about 50% of the entries.
    for (idx r = 0 ; r < testSize ; r++)
        for (idx c = 0 ; c < testSize ; c++){
            if (rand() % 2 ){
                real val = -1.0 + ( (real) (rand() % 1000) ) / 500.0; // RANDOM VALUE IN -1 -> +1 RANGE
                mm.addNonZero(r, c, val);
            }
        }
    SpaMtrix::IRCMatrix Atemp = mm.getIRCMatrix();

    // Assignment test.
    SpaMtrix::IRCMatrix A;
    A = Atemp;
    // Display the matrix and its sparsity pattern.
    A.print();
    A.spy();

    // Matrix-vector multiplication.
    std::cout << "Performing Ax = b, where x is all ones" << std::endl;
    SpaMtrix::Vector x(testSize);
    x = 1.0;
    SpaMtrix::Vector b = A * x;
    b.print("b");

    // Matrix-scalar operations.
    SpaMtrix::IRCMatrix M = -10*A;
    M*=-0.1;
    M.print();



    return 0;
}

