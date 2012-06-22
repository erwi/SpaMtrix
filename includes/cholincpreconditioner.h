#ifndef CHOLINCPRECONDITIONER_H
#define CHOLINCPRECONDITIONER_H
#include <setup.h>
#include <preconditioner.h>
#include <ircmatrix.h>
#include <vector.h>


class CholIncPreconditioner:Preconditioner
{

public:
  CholIncPreconditioner(const IRCMatrix& A);
  void solveMxb(Vector &x, const Vector &b) const;
};

#endif
