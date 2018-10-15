#CXX   := /usr/local/bin/g++-7

CXX   := g++
TAR   := tar xzf
MKDIR := mkdir -p
CD    := cd
CP    := cp
CMAKE := cmake

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

FILE  := lapack-3.7.1/make.inc
EXIST := $(shell ls | grep ${FILE})

.PHONY: all
all: mmseqs lapack gs2 clean

.PHONY: mmseqs
mmseqs:
	$(TAR) MMseqs2.tar.gz
	$(MKDIR) MMseqs2/build
	$(CD) MMseqs2/build; $(CMAKE) -DHAVE_SSE4_1=1 -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..
	$(MAKE) -C MMseqs2/build
	$(MAKE) -C MMseqs2/build install 

.PHONY: lapack
lapack:
	$(TAR) lapack-3.7.1.tar.gz 
	$(MKDIR) lib
    ifneq (${EXIST}, ${FILE})
	$(CP) lapack-3.7.1/make.inc.example lapack-3.7.1/make.inc
    endif
	$(MAKE) -C lapack-3.7.1 blaslib
	$(MAKE) -C lapack-3.7.1 cblaslib
	$(MAKE) -C lapack-3.7.1 lapacklib
	$(MAKE) -C lapack-3.7.1 lapackelib
	$(CP) -f lapack-3.7.1/*.a lib
	$(CP) -f lapack-3.7.1/CBLAS/include/*.h lib
	$(CP) -f lapack-3.7.1/LAPACKE/include/*.h lib

gs2: eigen.o transitivity.o $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIB)

eigen.o: eigen.cpp
	$(CXX) $(CXXFLAGS) $(INC) -c $<

transitivity.o: transitivity.cpp
	$(CXX) $(CXXFLAGS) $(INC) -c $<

$(OBJECTS): %.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<

.PHONY: clean
clean:
	$(RM) eigen.o transitivity.o $(OBJECTS)
