
include ../../makefile.def


pcg_test:pcg_test.o $(LIBRARIES_PATH)/libspamtrix.a
	g++ -o pcg_test pcg_test.o -L$(LIBRARIES_PATH) -lspamtrix -lgomp
	


pcg_test.o:pcg_test.cpp
	g++ -c pcg_test.cpp -I$(INCLUDE_PATH)
