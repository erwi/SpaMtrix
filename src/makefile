include ../makefile.def

OBJS = spamtrix_blas.o ircmatrix.o \
	vector.o\
	diagpreconditioner.o\
	cholincpreconditioner.o\
	matrixmaker.o  \
	tdmatrix.o\
	cholesky.o \
	lu.o\
	fleximatrix.o\
	writer.o\
	densematrix.o\
	iterativesolvers.o

LIBS += -lgomp


all:$(OBJS)
.cpp.o:
	g++ -c -g $< $(CXXFLAGS) -I$(INCLUDE_PATH) 

clean:
	rm -f *.o; rm -f *.exe
	cd ../lib; rm -f *.*

lib:$(OBJS)
	ar rcs libspamtrix.a $(OBJS)
	mv libspamtrix.a ../lib
