#ifndef LUINCPRECONDITIONER_H
#define LUINCPRECONDITIONER_H

#include <iostream>
#include <math.h>

#include <spamtrix_preconditioner.hpp>
#include <spamtrix_fleximatrix.hpp>
#include <spamtrix_ircmatrix.hpp>
#include <spamtrix_vector.hpp>

namespace SpaMtrix
{

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
} // end namespace SpaMtrix

#endif // LUINCPRECONDITIONER_H
