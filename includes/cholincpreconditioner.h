#ifndef CHOLINCPRECONDITIONER_H
#define CHOLINCPRECONDITIONER_H
#include <setup.h>
#include <ircmatrix.h>
#include <preconditioner.h>
#include <vector.h>

#include <iostream>
using std::cout;
using std::endl;

// THIS FILE IS ICLUDED FROM 'ircmatrix.h', SO NEED TO DECLARE IRCMatrix HERE
class IRCMatrix; 

class CholIncPreconditioner:public Preconditioner
{
  IRCMatrix* LD; // LOWER DIAGONAL INCOMPLETE CHOLESKY FACTORISED MATRIX

  void makeLowerDiagonal(const IRCMatrix& A);

public:
  CholIncPreconditioner(const IRCMatrix &A);
  ~CholIncPreconditioner() {if (LD) delete LD;}
  void solveMxb(Vector &x, const Vector &b) const;
};

#endif
