#ifndef PRECONDITIONER_H
#define PRECONDITIONER_H
#include <vector.h>


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
  
};

} // end namespace SpaMtrix

#endif
