include ../../makefile.def


tridiag_solver:tridiag_solver.o $(LIBRARIES_PATH)/libspamtrix.a
	g++ -o tridiag_solver tridiag_solver.o -L$(LIBRARIES_PATH) -lspamtrix -lgomp
	


tridiag_solver.o:tridiag_solver.cpp
	g++ -c tridiag_solver.cpp -I$(INCLUDE_PATH) 
