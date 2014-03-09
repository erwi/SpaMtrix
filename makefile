
clean:
	cd src; make clean
	cd lib; rm -f *.*


	cd examples/decomposition_solvers; make clean
	cd examples/pcg_test; make clean
	cd examples/test4; make clean
	cd examples/test5; make clean
	cd examples/test6; make clean
	cd examples/matrix_basics; make clean
	cd examples/input_output; make clean
	cd examples/tridiag_solver; make clean
	cd examples/eigenvalues; make clean
all:
	make libs
	make test


libs:
	@echo "LIBRARIES"
	cd src; make lib

python:
	@echo "MAKING PYTHIN WRAPPERS"
	cd src; make lib; make python


test:
	@echo "-------------------"
	@echo "  MATRIX BASICS    "
	@echo "-------------------"
	cd examples/matrix_basics; make
	@echo "-------------------"
	@echo "  INPUT OUTPUT     "
	@echo "-------------------"
	cd examples/input_output; make
	@echo "-------------------"
	@echo "     PCG TEST      "
	@echo "-------------------"
	cd examples/pcg_test; make
	@echo "-------------------"
	@echo "  TRIADIAG SOLVER  "
	@echo "-------------------"
	cd examples/tridiag_solver; make
	@echo "-----------------------"
	@echo "DECOMPOSITION SOLVERS :"
	@echo "-----------------------"
	cd examples/decomposition_solvers; make
	@echo "-------------------"
	@echo "      TEST 4:"
	@echo "-------------------"
	cd examples/test4; make test
	@echo "-------------------"
	@echo "      TEST 5:"
	@echo "-------------------"
	cd examples/test5; make test
	@echo "-------------------"
	@echo "      TEST 6:"
	@echo "-------------------"
	cd examples/test6; make test
	@echo "-------------------"
	@echo "   EIGENVALUES:"
	@echo "-------------------"
	cd examples/eigenvalues; make eigenvalues


