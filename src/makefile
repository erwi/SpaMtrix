
OBJS = spamtrix_blas.o ircmatrix.o ircmatrix_test.o makeircmatrix.o vector.o spamtrix_blas.o cg.o pcg.o diagpreconditioner.o

OBJS += $(IRCMATOBJS)
EXECS = ircmatrix_test.exe

CXXFLAGS = -Wall -DDEBUG

INCDIR = ../includes
all:$(OBJS)
.cpp.o:
	g++ -c -g $< -I$(INCDIR) $(CXXFLAGS)
	

ircmatrix_test.exe: $(IRCMATOBJS) $(OBJS)
	g++ -o $@ $(IRCMATOBJS)
	
clean:
	rm -f *.o; rm -f *.exe
	
	