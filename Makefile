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
all: mmseqs lapack gs2 clean

.PHONY: mmseqs
mmseqs:
	tar xzf MMseqs2.tar.gz
	mkdir -p MMseqs2/build
	cd MMseqs2/build; cmake -DHAVE_SSE4_1=1 -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..
	make -C MMseqs2/build
	make -C MMseqs2/build install 

.PHONY: lapack
lapack:
	tar xvzf lapack-3.7.1.tar.gz 
	mkdir -p lib
	cp lapack-3.7.1/make.inc.example lapack-3.7.1/make.inc
	make -C lapack-3.7.1 blaslib
	make -C lapack-3.7.1 cblaslib
	make -C lapack-3.7.1 lapacklib
	make -C lapack-3.7.1 lapackelib
	cp lapack-3.7.1/*.a lib
	cp lapack-3.7.1/CBLAS/include/*.h lib
	cp lapack-3.7.1/LAPACKE/include/*.h lib

gs2: eigen.o $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIB)

eigen.o: eigen.cpp
	$(CXX) $(CXXFLAGS) $(INC) -c $<

$(OBJECTS): %.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<

.PHONY: clean
clean:
	$(RM) eigen.o $(OBJECTS)
