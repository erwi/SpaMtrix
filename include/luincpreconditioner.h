#ifndef LUINCPRECONDITIONER_H
#define LUINCPRECONDITIONER_H

#include <preconditioner.h>
#include <fleximatrix.h>
#include <ircmatrix.h>
#include <vector.h>

#include <iostream>
#include <math.h>

class LUIncPreconditioner: public Preconditioner
{


    FlexiMatrix a;  // stores both L and U. L is unit lower diagonal, so its
                    // diagonal values are not stored

    LUIncPreconditioner(){}


    void forwardSubstitution(Vector &x, const Vector &b) const;
    void backwardSubstitution(Vector &x, const Vector &b) const;

public:
    LUIncPreconditioner(const IRCMatrix &A);
    void print() const;
    void solveMxb(Vector &x, const Vector &b) const;
};


#endif // LUINCPRECONDITIONER_H
