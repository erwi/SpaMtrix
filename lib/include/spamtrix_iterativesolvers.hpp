#ifndef ITERATIVESOLVERS_H
#define ITERATIVESOLVERS_H

#include <spamtrix_setup.hpp>
#include <spamtrix_ircmatrix.hpp>
#include <spamtrix_vector.hpp>
#include <spamtrix_preconditioner.hpp>
#include <spamtrix_densematrix.hpp>
namespace SpaMtrix
{

/**
 * @brief Iterative linear-system solvers and their control parameters.
 */
class IterativeSolvers
{
  
public:
    /** @brief Maximum iteration count. */
    idx maxIter;
    /** @brief Maximum inner iteration count used by GMRES. */
    idx maxInnerIter;
    /** @brief Requested or achieved numerical tolerance. */
	real toler;

    /** @brief Construct with default iteration limits. */
    IterativeSolvers();
    /**
     * @brief Construct with a maximum iteration count and tolerance.
     *
     * @param maxIter Maximum iterations.
     * @param toler Requested tolerance.
     */
    IterativeSolvers(const idx maxIter, const real toler);
    /**
     * @brief Construct with outer and inner iteration limits.
     *
     * @param maxIter Maximum outer iterations.
     * @param maxInnerIter Maximum inner iterations.
     * @param toler Requested tolerance.
     */
    IterativeSolvers(const idx maxIter, 
		     const idx maxInnerIter,
		     const real toler);
    
    
    /**
     * @brief Solve a linear system with preconditioned conjugate gradients.
     *
     * @param A System matrix.
     * @param x Initial guess on input; solution on output.
     * @param b Right-hand-side vector.
     * @param M Preconditioner.
     * @return True if convergence was achieved.
     */
    bool pcg( const IRCMatrix &A,
              Vector &x,
              const Vector &b,
              const Preconditioner &M
              );

    /**
     * @brief Solve a linear system with GMRES.
     *
     * @param A System matrix.
     * @param x Initial guess on input; solution on output.
     * @param b Right-hand-side vector.
     * @param M Preconditioner.
     * @return True if convergence was achieved.
     */
    bool gmres(const IRCMatrix &A,
                      Vector &x,
                      const Vector &b,
                      const Preconditioner &M);

};
} // end namespace SpaMtrix
#endif // ITERATIVESOLVERS_H
