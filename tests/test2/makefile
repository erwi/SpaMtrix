include ../../makefile.def


test:test2.o $(LIBRARIES_PATH)/libspamtrix.a
	g++ -o test2 test2.o -L$(LIBRARIES_PATH) -lspamtrix -lgomp
	


test2.o:test2.cpp
	g++ -c test2.cpp -I$(INCLUDE_PATH) 
