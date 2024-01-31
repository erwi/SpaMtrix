#ifndef CHOLINCPRECONDITIONER_H
#define CHOLINCPRECONDITIONER_H
#include <iostream>
#include <math.h>

#include <spamtrix_setup.hpp>
#include <spamtrix_ircmatrix.hpp>
#include <spamtrix_vector.hpp>
#include <spamtrix_fleximatrix.hpp>
#include <spamtrix_preconditioner.hpp>

namespace SpaMtrix {
class CholIncPreconditioner: public Preconditioner {
    FlexiMatrix L;
    CholIncPreconditioner(): L() {}
    void forwardSubstitution(Vector&x, const Vector& b) const;
    void backwardSubstitution(Vector&x, const Vector& b) const;
public:
    CholIncPreconditioner(const IRCMatrix &A);
    virtual ~CholIncPreconditioner();
    void print() const;
    void solveMxb(Vector &x, const Vector &b) const;
};
} // end namespace SpaMtrix
#endif
