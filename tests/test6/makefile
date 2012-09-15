include ../../makefile.def


test:test6.o $(LIBRARIES_PATH)/libspamtrix.a
	g++ -o test6 test6.o -L$(LIBRARIES_PATH) -lspamtrix -lgomp
	


test6.o:test6.cpp
	g++ -c test6.cpp -I$(INCLUDE_PATH) 
