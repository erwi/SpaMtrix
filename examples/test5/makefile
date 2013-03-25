include ../../makefile.def


test:test5.o $(LIBRARIES_PATH)/libspamtrix.a
	g++ -o test5 test5.o -L$(LIBRARIES_PATH) -lspamtrix -lgomp
	


test5.o:test5.cpp
	g++ -c test5.cpp -I$(INCLUDE_PATH) 

clean:
	rm -f *.csv
	rm -f *.o
	rm -f ./test5