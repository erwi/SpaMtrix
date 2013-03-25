/*!
 * Test file for Cholesky and LU decompositions. A 5-point Poisson finite differences
 * problem is solved using both Cholesky and LU decomposition solvers.
 * 
 * Use: ./decomposition_solvers gridLen numThreads,\n
 * where gridLen specifies FD grid side length and \n
 * numThreads the number of threads used by OpenMP.\n
 * 
 */


#include <iostream>
#include <chrono>
#include <omp.h>

#include <spamtrix_setup.hpp>
#include <spamtrix_matrixmaker.hpp>
#include <spamtrix_ircmatrix.hpp>
#include <spamtrix_lu.hpp>
#include <spamtrix_cholesky.hpp>
#include <spamtrix_blas.hpp>
#include <spamtrix_tickcounter.hpp>


using namespace SpaMtrix;
using std::cout;
using std::endl;
int main(int nargs, char* args[])
{
    if (nargs == 1){
        cout << "arguments: \nfirst = grid size (default is 5)" << endl;
        cout << "second = number of threads (default is 1)" << endl;
    }
    
    idx gridLen = 5;
    idx numThreads = 1;
    if (nargs>1)
        gridLen = atoi(args[1]);

    if (nargs>2)
        numThreads = atoi(args[2]);

    idx numDoF = gridLen*gridLen;
    MatrixMaker mm(numDoF,numDoF);
    mm.poisson5Point();

    IRCMatrix A = mm.getIRCMatrix();

    cout << "\nMatrix size is : " << numDoF << "x" << numDoF << endl;
    // PERFORMANCE TIMER
    TickCounter<std::chrono::milliseconds> timer;

    // SET NUMBER OF THREAD TO USE
    omp_set_num_threads(numThreads);
    cout << "num theads : " << numThreads << endl <<endl;

    // CREATE LU DECOMPOSITION
    timer.start();
    LU lu(A);
    timer.stop();
    size_t t = timer.getElapsed();
    cout << "LU decomposition time [ms] : " << t << endl;
    timer.reset();

    // CREATE CHOLESKY DECOMPOSITION
    timer.start();
    Cholesky C(A);
    timer.stop();
    size_t tch = timer.getElapsed();
    timer.reset();
    cout << "Cholesky decomposition time [ms] : " << tch << endl<< endl;

    
    // SOLUTION VECTORS X FOR LU AND CHOLESKY
    Vector xlu(numDoF);
    Vector xch(numDoF);
    Vector b(numDoF);

    // RHS, Ax = 1
    b.setAllValuesTo(1.0);

    // SOLVE LU
    timer.start();
    lu.solve(xlu,b);
    timer.stop();
    size_t tlus = timer.getElapsed();
    timer.reset();
    cout << "LU solution time [ms] : " << tlus << endl;


    // SOLVE CHOLESKY
    timer.start();
    C.solve(xch,b);
    timer.stop();
    size_t tcs = timer.getElapsed();
    cout << "Cholseky solution time [ms] : " << tcs << endl;


    // MAKE SURE BOTH METHODS GIVE SAME(ISH) ANSWER
    real diff = norm(xlu-xch);
    cout << "solution difference : " << diff << endl;
    return 0;
}

