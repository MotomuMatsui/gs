CXX := g++

VPATH := src
INC   := -Ilib
LIB   := -Llib -llapacke -llapack -lcblas -lrefblas -lgfortran -lm

CXXFLAGS := -O3
CXXFLAGS += -std=c++11
CXXFLAGS += -march=native
CXXFLAGS += -fno-exceptions
CXXFLAGS += -Wall

OBJECTS  := messages.o
OBJECTS  += format.o
OBJECTS  += mmseqs.o
OBJECTS  += gs_functions.o
OBJECTS  += sc_functions.o
OBJECTS  += sc.o
OBJECTS  += ep.o
OBJECTS  += gs.o
OBJECTS  += main.o

.PHONY: all
all: gs2 clean

gs2: eigen.o $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(INC) -o $@ $^ $(LIB)	

eigen.o: eigen.cpp
	$(CXX) $(CXXFLAGS) $(INC) -c $< $(LIB)

$(OBJECTS): %.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<

.PHONY: clean
clean:
	$(RM) eigen.o $(OBJECTS)
