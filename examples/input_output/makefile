
include ../../makefile.def


input_output:input_output.o $(LIBRARIES_PATH)/libspamtrix.a
	g++ -o input_output input_output.o -L$(LIBRARIES_PATH) -lspamtrix -lgomp
	


input_output.o:input_output.cpp
	g++ -c input_output.cpp -I$(INCLUDE_PATH)

clean:
	rm -f *.o
	rm -f ./input_output
	rm -f *.txt