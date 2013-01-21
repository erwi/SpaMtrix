#ifndef DIAGPRECONDITIONER_H
#define DIAGPRECONDITIONER_H
#include <omp.h>

#include <spamtrix_setup.hpp>
#include <spamtrix_vector.hpp>
#include <spamtrix_ircmatrix.hpp>
#include <spamtrix_preconditioner.hpp>


namespace SpaMtrix
{
class DiagPreconditioner: public Preconditioner {
  
    Vector diagonal;
    
    DiagPreconditioner();	// PRIVATE COSTRUCTOR
public:
  DiagPreconditioner(const IRCMatrix& A);
   
  //void applyToVector(Vector& v) const;
  void solveMxb(Vector &x, const Vector &b) const;  
  
};
} // end namespace SpaMtrix

#endif

