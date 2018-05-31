CXX := g++
INC := -Ilib
LIB := -Llib -llapacke -llapack -lcblas -lrefblas -lgfortran -lm
CXXFLAGS += -O3
CXXFLAGS += -std=c++11
CXXFLAGS += -march=native
CXXFLAGS += -fno-exceptions
CXXFLAGS += -Wall

OBJECTS = messages.o format.o eigen.o mmseqs.o gs_functions.o sc_functions.o sc.o ep.o gs.o main.o

all: gs2

gs2: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(INC) -o $@ $^ $(LIB)	

$(OBJECTS): %.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INC) -c $< $(LIB)

clean:
	$(RM) gs2
	$(RM) $(OBJECTS)