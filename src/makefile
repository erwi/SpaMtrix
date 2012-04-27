
OBJS = spamtrix_blas.o ircmatrix.o vector.o cg.o pcg.o diagpreconditioner.o matrixmaker.o test1.o test2.o

OBJS += $(IRCMATOBJS)
EXECS = ircmatrix_test.exe

CXXFLAGS += -Wall 
CXXFLAGS += -DUSES_OPENMP # WHETHER COMPILED WITH OPENMP FOR LOOPS
#CXXFLAGS += -DDEBUG -g
CXXFLAGS += -O3 -g -fopenmp -march=native -ffast-math

INCDIR = ../includes

LIBS += -lgomp


all:$(OBJS)
.cpp.o:
	g++ -c -g $< -I$(INCDIR) $(CXXFLAGS)

clean:
	rm -f *.o; rm -f *.exe
	

TEST1_OBJS = test1.o ircmatrix.o vector.o cg.o spamtrix_blas.o matrixmaker.o
test1:
	g++ -o test1.exe $(TEST1_OBJS) $(LIBS)

TEST2_OBJS = test2.o ircmatrix.o vector.o diagpreconditioner.o cg.o pcg.o spamtrix_blas.o matrixmaker.o	
test2:
	g++ -o test2.exe $(TEST2_OBJS) $(LIBS)