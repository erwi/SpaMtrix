#include <diagpreconditioner.h>

DiagPreconditioner::DiagPreconditioner(const IRCMatrix& A ):
diagonal(A.getNumRows() )
{
    for (idx i = 0 ; i < diagonal.getLength() ; ++i)
    {
      diagonal[i] = A.sparse_get(i,i);
    }
}

/*
void DiagPreconditioner::applyToVector(Vector& v) const
{
/*!
 * Performs v = M^(-1)v
 */  

//#ifdef DEBUG
//  assert(diagonal.getLength() == v.getLength() );
//#endif
//#pragma omp parallel for
//    for (idx i = 0 ; i < diagonal.getLength() ; ++i)
//      v[i]/=diagonal[i];
//}
  
void DiagPreconditioner::solveMxb(Vector &x, Vector &b) const
{
#ifdef DEBUG
  assert(diagonal.getLength() == x.getLength() );
  assert(diagonal.getLength() == b.getLength() );
#endif
  
#ifdef USES_OPENMP
#pragma omp parallel for
#endif
  for (idx i = 0 ; i < diagonal.getLength() ; ++i)
    x[i] = b[i] / diagonal[i];
  
}