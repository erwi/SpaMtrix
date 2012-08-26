include ../makefile.def

OBJS = spamtrix_blas.o ircmatrix.o \
	vector.o cg.o pcg.o\
	diagpreconditioner.o\
	cholincpreconditioner.o\
	matrixmaker.o  tdmatrix.o\
	cholesky.o fleximatrix.o\
	writer.o

LIBS += -lgomp


all:$(OBJS)
.cpp.o:
	g++ -c -g $< -I$(INCLUDE_PATH) $(CXXFLAGS)

clean:
	rm -f *.o; rm -f *.exe
	cd ../lib; rm -f *.*

lib:$(OBJS)
	ar rcs libspamtrix.a $(OBJS)
	mv libspamtrix.a ../lib

#TEST1_OBJS = test1.o ircmatrix.o vector.o cg.o #spamtrix_blas.o matrixmaker.o
#test1:
#	g++ -o test1.exe $(TEST1_OBJS) $(LIBS)

#TEST2_OBJS = test2.o ircmatrix.o vector.o #diagpreconditioner.o cg.o pcg.o spamtrix_blas.o #matrixmaker.o tdmatrix.o
#test2:
#	g++ -o test2.exe $(TEST2_OBJS) $(LIBS)

#TEST3_OBJS = test3.o ircmatrix.o matrixmaker.o cholincpreconditioner.o
#test3:
#	g++ -o test3.exe $(TEST3_OBJS) $(LIBS)