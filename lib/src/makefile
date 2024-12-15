include ../makefile.def

OBJS = spamtrix_blas.o\
	ircmatrix.o \
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
	reader.o\
	densematrix.o\
	iterativesolvers.o\
	powermethod.o

LIBS += -lgomp

all:$(OBJS)
.cpp.o:
	g++ -c -g $< $(CXXFLAGS) -I$(INCLUDE_PATH) 

clean:
	rm -f *.o; rm -f *.exe
	rm -f *.py; rm -f *.pyc; rm -f spamtrix_wrap.*; rm -f *.so
	cd ../lib; rm -f *.*

lib:$(OBJS)
	ar rcs libspamtrix.a $(OBJS)
	mkdir -p ../lib
	mv libspamtrix.a ../lib

# USE SWIG TO MAKE PYTHON BINDINGS
python:
	swig -c++ -v -python spamtrix.i
	g++ -fPIC -c spamtrix_wrap.cxx -I/usr/include/python2.7 -I$(INCLUDE_PATH)
	g++ -shared $(OBJS) spamtrix_wrap.o -o _SpaMtrix.so
	mv _SpaMtrix.so ../lib
	mv SpaMtrix.py ../lib