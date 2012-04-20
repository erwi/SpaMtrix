
OBJS = mvvd.o rcmatrix.o 
IRCMATOBJS = ircmatrix.o ircmatrix_test.o makeircmatrix.o
OBJS += $(IRCMATOBJS)
EXECS = ircmatrix_test.exe

CXXFLAGS = -Wall -DDEBUG

INCDIR = ../includes
all:$(OBJS)
.cpp.o:
	g++ -c -g $< -I$(INCDIR) $(CXXFLAGS)
	

ircmatrix_test.exe: $(IRCMATOBJS)
	g++ -o $@ $(IRCMATOBJS)
	
clean:
	rm -f *.o; rm -f *.exe
	
	