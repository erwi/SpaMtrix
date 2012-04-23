
OBJS = spamtrix_blas.o ircmatrix.o test1.o vector.o cg.o pcg.o diagpreconditioner.o matrixmaker.o

OBJS += $(IRCMATOBJS)
EXECS = ircmatrix_test.exe

CXXFLAGS = -Wall -DDEBUG

INCDIR = ../includes
all:$(OBJS)
.cpp.o:
	g++ -c -g $< -I$(INCDIR) $(CXXFLAGS)

clean:
	rm -f *.o; rm -f *.exe
	

TEST1_OBJS = test1.o ircmatrix.o vector.o cg.o spamtrix_blas.o matrixmaker.o
test1:
	g++ -o test1.exe $(TEST1_OBJS)
	