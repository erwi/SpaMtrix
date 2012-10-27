include ../makefile.def

OBJS = spamtrix_blas.o ircmatrix.o \
	vector.o\
	diagpreconditioner.o\
	cholincpreconditioner.o\
	luincpreconditioner.o\
	matrixmaker.o  \
	tdmatrix.o\
	cholesky.o \
	lu.o\
	fleximatrix.o\
	writer.o\
	densematrix.o\
	iterativesolvers.o

LIBS += -lgomp


PYSPAMTRIX_OBJS = vector.o


all:$(OBJS)
.cpp.o:
	g++ -c -g $< $(CXXFLAGS) -I$(INCLUDE_PATH) 

clean:
	rm -f *.o; rm -f *.exe
	rm -f *.py; rm -f *.pyc; rm -f spamtrix_wrap.*; rm -f *.so
	cd ../lib; rm -f *.*

lib:$(OBJS)
	ar rcs libspamtrix.a $(OBJS)
	mv libspamtrix.a ../lib

# USE SWIG TO MAKE PYTHON BINDINGS
pyspamtrix:
	swig -c++ -v -python spamtrix.i
	g++ -fPIC -c spamtrix_wrap.cxx -I/usr/include/python2.7 -I$(INCLUDE_PATH)
	g++ -fPIC -c vector.cpp $(CXXFLAGS) -I$(INCLUDE_PATH) -I/usr/include/python2.7
	g++ -shared $(PYSPAMTRIX_OBJS) spamtrix_wrap.o -o _SpaMtrix.so