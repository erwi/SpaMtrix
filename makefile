 
clean:
	cd src; make clean
	cd lib; rm -f *.*
	cd examples/test1; rm -f *.o
	cd examples/test2; rm -f *.o
	cd examples/test3; rm -f *.o
	cd examples/test4; rm -f *.o
	cd examples/test5; rm -f *.o
	cd examples/test6; rm -f *.o
	cd examples/matrix_basics; make clean 
	cd examples/input_output; make clean
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
	@echo "      TEST 2:"
	@echo "-------------------"
	cd examples/test2; make test
	@echo "-------------------"
	@echo "      TEST 3:"
	@echo "-------------------"
	cd examples/test3; make test
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


