#ifndef DIAGPRECONDITIONER_H
#define DIAGPRECONDITIONER_H
#include <setup.h>
#include <vector.h>
#include <ircmatrix.h>
#include <preconditioner.h>
#include <omp.h>

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

