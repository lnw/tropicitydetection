
CXX=g++

FLAGS=-O2 -Wshadow -Wunused -Wpedantic -std=c++11 -fPIC
PYTHON_INCLUDES_C=-I/usr/include/python3.7m

HEADERS=geometry3.hh auxiliary.hh trajectory.hh cube.hh
OBJECTS=geometry3.o trajectory.o cube.o
OBJECTS_STANDALONE=main.o
OBJECTS_LIB=interface.o
OBJECTS_P=$(patsubst %.o, build/%.o, $(OBJECTS))
OBJECTS_STANDALONE_P=$(patsubst %.o, build/%.o, $(OBJECTS_STANDALONE))
OBJECTS_LIB_P=$(patsubst %.o, build/%.o, $(OBJECTS_LIB))

build/%.o: %.cc $(HEADERS) Makefile
	$(CXX) $(FLAGS) $(PYTHON_INCLUDES_C) -c $< -o $@

main: $(OBJECTS_P) main.cc Makefile
	$(CXX) $(FLAGS) $(OBJECTS_P) main.cc -o $@

libtropicity.so: $(OBJECTS_P) $(OBJECTS_LIB_P) Makefile
	$(CXX) $(FLAGS) -shared $(OBJECTS_P) $(OBJECTS_LIB_P) -o $@



.PHONY: clean
clean:
	rm -f a.out build/*.o
