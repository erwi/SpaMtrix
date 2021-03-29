#ifndef POWERMETHOD_H
#define POWERMETHOD_H
#include <spamtrix_blas.hpp>
#include <spamtrix_setup.hpp>
namespace SpaMtrix{
idx powerMethod(const SpaMtrix::IRCMatrix &A,
                 real &eigenValue,
                 SpaMtrix::Vector &eigenVector,
                 real &toler,
                 unsigned int maxIter = MAX_INDEX
                );

//{
/*!
 * Calculates the dominant eigenvalue and corresponding eigenvector of matrix A
 * using power iterations. \n
 * The eigenvalue and aigenvector are returned in parameters 'eigenValue' and
 * 'eigenVector', respectively.
 * 'toler' specifies the required accuracy of the method, as a percentage of the 
 * calculated eigenvalue (how much the eigenvalue changes during each iteration relative to its 
 * value). \n
 * The function terminates if more than 'maxIter' iterations are used, and the total number of iterations used is returned. 
 */
/*
    // MAKE SURE INITIAL GUESS IS NOT A ZERO VECTOR
    real n = eigenVector.getNorm();
    if (n == 0.0){
        eigenVector(0) = 1.0;
    }
    eigenVector.normalise();
    
    real dL = toler+1.0; // DELTA EIGENVALUE
    //real Lo(0);          // PREVIOUS EIGENVALUE
        
    SpaMtrix::Vector q(eigenVector);
    SpaMtrix::Vector z = A*eigenVector;
    idx iter(0);
    // DO WHILE RELATIVE CHANGE IS LESS THAN toler
    while ( dL > toler){
       z.normalise();
       q = z;
       z = A*q;
       
       // UPDATE EIGENVALUE AND CALCULATE RELATIVE CHANGE
       real Lo = eigenValue;    // OLD
       eigenValue = dot(q,z);   // eigenValue = xAx
       dL = fabs(Lo-eigenValue) / eigenValue;
       
       iter++;
       if (iter >= maxIter ){
           break;
       }
    }
    eigenVector = q;
    return iter;
}
*/
} // end namespace SpaMtrix
#endif
