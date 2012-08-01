 
clean:
	cd src; make clean
	cd tests/test1; rm *.o
	cd tests/test2; rm *.o
	cd tests/test4; rm *.o

all:
	make libs
	make test


libs:
	@echo "LIBRARIES"
	cd src; make lib

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
	@echo "      TEST 4:"
	@echo "-------------------"
	cd tests/test4; make test