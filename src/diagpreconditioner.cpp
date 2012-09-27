#include <diagpreconditioner.h>

DiagPreconditioner::DiagPreconditioner(const IRCMatrix& A ):
diagonal(A.getNumRows() )
{
    for (idx i = 0 ; i < diagonal.getLength() ; ++i)
    {
      diagonal[i] = A.sparse_get(i,i);
    }
}

  
void DiagPreconditioner::solveMxb(Vector &x, const Vector &b) const
{
#ifdef DEBUG
  assert(diagonal.getLength() == x.getLength() );
  assert(diagonal.getLength() == b.getLength() );
#endif
  
#ifdef USES_OPENMP
//#pragma omp parallel for schedule(static,1000)
#endif
  for (idx i = 0 ; i < diagonal.getLength() ; ++i)
    x[i] = b[i] / diagonal[i];
  
}