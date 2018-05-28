LAPACK_VER=3.8.0
LIB_PATH=../lib

# Decompress LAPACK package
tar xvzf lapack-${LAPACK_VER}.tar.gz
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

# copy Library Files
cp *.a ${LIB_PATH}

# copy Header Files
cp CBLAS/include/*.h ${LIB_PATH}
cp LAPACKE/include/*.h ${LIB_PATH}

cd ..