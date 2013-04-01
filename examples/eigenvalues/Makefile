include ../../makefile.def 
#SPATH = /home/eero/Projects/SpaMtrix/
#LPATH = $(SPATH)lib
#IPATH = $(SPATH)include
eigenvalues:eigenvalues.o
	g++ -o eigenvalues eigenvalues.o -L$(LIBRARIES_PATH) -lspamtrix

eigenvalues.o:eigenvalues.cpp
	g++ -c eigenvalues.cpp  -std=c++0x -I$(INCLUDE_PATH)

clean:
	rm -f *.o
	rm -f ./eigenvalues
	rm -f ./*.csv