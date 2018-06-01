#!/usr/bin/sh

LAPACK_VER=3.8.0
#LAPACK_VER=3.7.1
LIB_PATH=../lib

# Decompress LAPACK package
tar xvzf lapack-${LAPACK_VER}.tar.gz
#tar xvzf lapack-${LAPACK_VER}.tgz
if [ ! -d lib ]; then mkdir lib; fi
cd lapack-${LAPACK_VER}

# Copy Makefile
cp make.inc.example make.inc

# Compile BLAS
make blaslib

# Compile CBLAS
make cblaslib

# Compile LAPACK
make lapacklib

# Compile LAPACKE
make lapackelib

# Copy Library Files
cp *.a ${LIB_PATH}

# Copy Header Files
cp CBLAS/include/*.h ${LIB_PATH}
cp LAPACKE/include/*.h ${LIB_PATH}

cd ..