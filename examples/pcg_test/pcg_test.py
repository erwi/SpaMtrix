"""@package docstring

 Python test\n
 
 Example/test showing how to create a sparse matrix using the MatrixMaker.\n
 
 A small 2x2 matrix problem Ax=b is constructed and solver using the 
 preconditioned conjugate gradient method.
  
 A = [[3,2];[2,6]] , b = [2,-8] is used.
"""
import SpaMtrix

mm = SpaMtrix.MatrixMaker(2,2)
mm.addNonZero(0,0,3) 
mm.addNonZero(0,1,2)
mm.addNonZero(1,0,2) 
mm.addNonZero(1,1,6)

A = SpaMtrix.IRCMatrix()
mm.makeSparseMatrix(A)
print "Solving Ax=b, where\nA is:"
A._print()

# RIGHT HAND SIDE VECTOR
b = SpaMtrix.Vector(2)
b[0]=2
b[1]=-8
print "b is : "
b._print("b")
# UNKNOWN VECTOR x
x = SpaMtrix.Vector(2)


# CRATE DIAGONAL PRECONDITIONER MATRIX
M = SpaMtrix.DiagPreconditioner(A)

# ITERATIVE SOLUTION
isol = SpaMtrix.IterativeSolvers(20,1e-7)
conv = isol.pcg(A,x,b,M)

strout = "convergence : "
if conv:
   strout +="OK"
else:
   strout +="NO"
   
print strout
print "maxIter = ", isol.maxIter
print "toler = " , isol.toler
print "solution vector : "
x._print("x")