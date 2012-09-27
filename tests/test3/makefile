include ../../makefile.def


test:test3.o $(LIBRARIES_PATH)/libspamtrix.a
	g++ -o test3 test3.o -L$(LIBRARIES_PATH) -lspamtrix -lgomp
	


test3.o:test3.cpp
	g++ -c test3.cpp -I$(INCLUDE_PATH) $(CXXFLAGS)
