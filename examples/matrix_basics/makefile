
include ../../makefile.def


matrix_basics:matrix_basics.o $(LIBRARIES_PATH)/libspamtrix.a
	g++ -o matrix_basics matrix_basics.o -L$(LIBRARIES_PATH) -lspamtrix -lgomp
	


matrix_basics.o:matrix_basics.cpp
	g++ -c matrix_basics.cpp -I$(INCLUDE_PATH)

clean:
	rm -f *.o
	rm -f ./matrix_basics