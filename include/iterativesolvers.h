#ifndef ITERATIVESOLVERS_H
#define ITERATIVESOLVERS_H

#include <setup.h>
#include <ircmatrix.h>
#include <vector.h>
#include <preconditioner.h>

class IterativeSolvers
{

public:

    IterativeSolvers();


    static bool pcg( const IRCMatrix &A,
              Vector &x,
              const Vector &b,
              const Preconditioner &M,
              idx &maxIter,
              real &toler
              );


    static bool gmres(const IRCMatrix &A,
                      Vector &x,
                      const Vector &b,
                      const Preconditioner &M,
                      idx &maxIter,
                      const idx &maxInnerIter,
                      real &toler
                      );

};

#endif // ITERATIVESOLVERS_H
