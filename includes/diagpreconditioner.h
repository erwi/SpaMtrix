#ifndef DIAGPRECONDITIONER_H
#define DIAGPRECONDITIONER_H
#include <setup.h>
#include <vector.h>
#include <ircmatrix.h>
#include <preconditioner.h>
#include <omp.h>
class DiagPreconditioner: public Preconditioner {
  
    Vector diagonal;
    
    DiagPreconditioner();	// PRIVATE COSTRUCTOR
public:
  DiagPreconditioner(const IRCMatrix& A);
   
  //void applyToVector(Vector& v) const;
  virtual void solveMxb(Vector &x, Vector &b) const;  
  
};


#endif

