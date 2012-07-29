
include ../../makefile.def


test:test1.o $(LIBRARIES_PATH)/libspamtrix.a
	g++ -o test1 test1.o -L$(LIBRARIES_PATH) -lspamtrix -lgomp
	


test1.o:test1.cpp
	g++ -c test1.cpp -I$(INCLUDE_PATH)
	