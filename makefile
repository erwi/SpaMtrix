 
clean:
	cd src; make clean

libs:
	@echo "LIBRARIES"
	cd src; make lib

test:
	cd tests/test1; make test
	cd tests/test2; make test
	@echo "-------------------"
	@echo "      TEST 4:"
	@echo "-------------------"
	cd tests/test4; make test