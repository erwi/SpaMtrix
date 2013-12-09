#ifndef PRECONDITIONER_H
#define PRECONDITIONER_H
#include <spamtrix_vector.hpp>


namespace SpaMtrix
{

// PURE VIRTUAL PRECONDITIONER BASE CLASS
class Preconditioner
{
public:
    virtual void solveMxb(Vector &x, const Vector &b) const  = 0; // SOLVES Mx = b

    Vector solve(const Vector &b) const
    {
        Vector x(b.getLength() );
        solveMxb( x , b);
        return x;
    }

    virtual ~Preconditioner();
};

} // end namespace SpaMtrix

#endif
