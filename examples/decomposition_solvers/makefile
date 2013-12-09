include ../../makefile.def


decomposition_solvers:decomposition_solvers.o $(LIBRARIES_PATH)/libspamtrix.a
	g++ -o decomposition_solvers decomposition_solvers.o -L$(LIBRARIES_PATH) -lspamtrix -lgomp

decomposition_solvers.o:decomposition_solvers.cpp
	g++ -c decomposition_solvers.cpp -I$(INCLUDE_PATH) $(CXXFLAGS)
	
clean:
	rm -f ./decomposition_solvers
	rm -f *.o
	
