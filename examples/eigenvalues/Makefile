include ../../makefile.def 
#SPATH = /home/eero/Projects/SpaMtrix/
#LPATH = $(SPATH)lib
#IPATH = $(SPATH)include
eigenvalues:eigenvalues.cpp $(LIBRARIES_PATH)/libspamtrix.a
	g++ -o eigenvalues eigenvalues.cpp -I$(INCLUDE_PATH) \
	-L$(LPATHIBRARIES_PATH) -lspamtrix -lgomp -std=c++0x

clean:
	rm -f *.o
	rm -f ./eigenvalues
	rm -f ./*.csv