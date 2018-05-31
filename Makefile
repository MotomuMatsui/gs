CXX := g++
INC := -Ilib
LIB := -Llib -llapacke -llapack -lcblas -lrefblas -lgfortran -lm
CXXFLAGS += -O3
CXXFLAGS += -std=c++11
CXXFLAGS += -march=native
CXXFLAGS += -fno-exceptions
CXXFLAGS += -Wall

all: gs2

gs2: messages.o format.o eigen.o mmseqs.o gs_functions.o sc_functions.o sc.o ep.o gs.o main.o
	$(CXX) $(CXXFLAGS) $(INC) -o $@ $^ $(LIB)	

messages.o: messages.cpp
	$(CXX) $(CXXFLAGS) $(INC) -c $< $(LIB)

format.o: format.cpp
	$(CXX) $(CXXFLAGS) $(INC) -c $< $(LIB)

eigen.o: eigen.cpp
	$(CXX) $(CXXFLAGS) $(INC) -c $< $(LIB)

mmseqs.o: mmseqs.cpp
	$(CXX) $(CXXFLAGS) $(INC) -c $< $(LIB)

gs_functions.o: gs_functions.cpp
	$(CXX) $(CXXFLAGS) $(INC) -c $< $(LIB)

sc_functions.o: sc_functions.cpp
	$(CXX) $(CXXFLAGS) $(INC) -c $< $(LIB)

sc.o: sc.cpp
	$(CXX) $(CXXFLAGS) $(INC) -c $< $(LIB)

ep.o: ep.cpp
	$(CXX) $(CXXFLAGS) $(INC) -c $< $(LIB)

gs.o: gs.cpp
	$(CXX) $(CXXFLAGS) $(INC) -c $< $(LIB)

main.o: main.cpp
	$(CXX) $(CXXFLAGS) $(INC) -c $< $(LIB)

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
	$(RM) gs2
