#ifndef DIAGPRECONDITIONER_H
#define DIAGPRECONDITIONER_H
#include <setup.h>
#include <vector.h>
#include <ircmatrix.h>
class DiagPreconditioner{
  
    Vector diagonal;
    
    DiagPreconditioner();	// PRIVATE COSTRUCTOR
public:
  DiagPreconditioner(const IRCMatrix& A);
   
  void applyToVector(Vector& v) const;
    
  
};


#endif

