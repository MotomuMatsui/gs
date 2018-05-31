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

all: gs2 cleaner

gs2: eigen.o $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(INC) -o $@ $^ $(LIB)	

eigen.o: eigen.cpp
	$(CXX) $(CXXFLAGS) $(INC) -c $< $(LIB)

$(OBJECTS): %.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<

cleaner:
	$(RM) eigen.o
	$(RM) $(OBJECTS)

clean:
	$(RM) gs2
	$(RM) eigen.o
	$(RM) $(OBJECTS)
