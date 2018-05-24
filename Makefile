CXX := g++
INC := -I./lib
LIB := -L./lib -llapacke -llapack -lcblas -lrefblas -lgfortran -lm
#LIB := -L./lib -llapacke -llapack -lcblas -lrefblas -L/usr/lib64 -lgfortran -lm
#LIB := -L./lib -llapacke -llapack -lcblas -lrefblas -lm
CXXFLAGS += -O3
#CXXFLAGS += -std=c++1z
CXXFLAGS += -std=c++11
CXXFLAGS += -march=native
CXXFLAGS += -fno-exceptions
CXXFLAGS += -Wall

all: gs

gs: messages.cpp format.cpp eigen.cpp mmseqs.cpp gs_functions.cpp sc_functions.cpp sc.cpp ep.cpp gs.cpp main.cpp
	$(CXX) $(CXXFLAGS) $(INC) -o $@ $^ $(LIB)	

clean:
	$(RM) main
	$(RM) format
	$(RM) mmseqs
	$(RM) gs
	$(RM) ep
	$(RM) sc
	$(RM) gs_functions
	$(RM) sc_functions
	$(RM) eigen
	$(RM) messages
