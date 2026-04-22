 
#include <iostream>
#include <spamtrix_ircmatrix.hpp>
#include <spamtrix_vector.hpp>
#include <spamtrix_matrixmaker.hpp>
#include <spamtrix_powermethod.hpp>
#include <spamtrix_blas.hpp>
#include <spamtrix_writer.hpp>
using std::endl;
using std::cout;


int main(int nargs, char* args[]){
    
    unsigned int n = 5;
    
    if (nargs >= 2 ){
        n = atoi(args[1]);
    }
    
    
    unsigned int m = n*n;
    SpaMtrix::MatrixMaker mm(m,m);
    mm.poisson5Point();
    SpaMtrix::IRCMatrix A = mm.getIRCMatrix();
  
    // Initial guess vector.
    SpaMtrix::Vector x(m);
    x(0) = 1.0;
          
    
    // Compute the dominant eigenvalue and eigenvector using power iterations.
    double L(0.0);
    double toler = 1e-12;
    idx iters = powerMethod(A,L,x, toler);
    
    x = abs(x);
    
    cout << iters << " iterations used " << endl;
    cout << "eigenValue : "<< L << endl;
    // Print the eigenvector if the problem is small enough.
    if (m <= 25){
        cout << "eigenVector e: " << endl;
        x.print("e");
    }
    
    SpaMtrix::Writer::writeCSV("out.csv",x, n);
    cout << "OK" << endl;
    return 0;
}