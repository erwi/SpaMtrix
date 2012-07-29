include ../../makefile.def


test:test4.o $(LIBRARIES_PATH)/libspamtrix.a
	g++ -o test4 test4.o -L$(LIBRARIES_PATH) -lspamtrix -lgomp
	


test4.o:test4.cpp
	g++ -c test4.cpp -I$(INCLUDE_PATH) 
