 
clean:
	cd src; make clean
	cd lib; rm -f *.*
	cd tests/test1; rm -f *.o
	cd tests/test2; rm -f *.o
	cd tests/test3; rm -f *.o
	cd tests/test4; rm -f *.o
	cd tests/test5; rm -f *.o
	cd tests/test6; rm -f *.o

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
	@echo "      TEST 1:"
	@echo "-------------------"
	cd tests/test1; make test
	@echo "-------------------"
	@echo "      TEST 2:"
	@echo "-------------------"
	cd tests/test2; make test
	@echo "-------------------"
	@echo "      TEST 3:"
	@echo "-------------------"
	cd tests/test3; make test
	@echo "-------------------"
	@echo "      TEST 4:"
	@echo "-------------------"
	cd tests/test4; make test
	@echo "-------------------"
	@echo "      TEST 5:"
	@echo "-------------------"
	cd tests/test5; make test
	@echo "-------------------"
	@echo "      TEST 6:"
	@echo "-------------------"
	cd tests/test6; make test


