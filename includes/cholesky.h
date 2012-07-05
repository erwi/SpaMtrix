#ifndef CHOLESKY_H
#define CHOLESKY_H
#include <setup.h>
#include <ircmatrix.h>
#include <vector>
#include <math.h>
#include <fleximatrix.h>
//struct ColVal; // FORWARD DECLARATION

class Cholesky
{

  FlexiMatrix L;

  Cholesky():L(0,0){}

public:
  Cholesky(const IRCMatrix& A);
  void print()const;
};


#endif // CHOLESKY_H
