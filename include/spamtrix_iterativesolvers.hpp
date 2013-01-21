#ifndef ITERATIVESOLVERS_H
#define ITERATIVESOLVERS_H

#include <spamtrix_setup.hpp>
#include <spamtrix_ircmatrix.hpp>
#include <spamtrix_vector.hpp>
#include <spamtrix_preconditioner.hpp>
#include <spamtrix_densematrix.hpp>
namespace SpaMtrix
{

class IterativeSolvers
{
  
public:
    // following variables act as input/output parameters
    idx maxIter;	// maximum/needed iteration count
    idx maxInnerIter;   // maximim allowed inner iterations
    real toler;		// required/achieved numerical accuracy

    IterativeSolvers();
    IterativeSolvers(const idx maxIter, const real toler);
    IterativeSolvers(const idx maxIter, 
		     const idx maxInnerIter,
		     const real toler);
    
    
    bool pcg( const IRCMatrix &A,
              Vector &x,
              const Vector &b,
              const Preconditioner &M
              );

    bool gmres(const IRCMatrix &A,
                      Vector &x,
                      const Vector &b,
                      const Preconditioner &M);

};
} // end namespace SpaMtrix
#endif // ITERATIVESOLVERS_H
