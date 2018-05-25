LAPACK_VER=3.8.0
LIB_PATH=../lib

# Download and decompression LAPACK package.
wget http://www.netlib.org/lapack/lapack-${LAPACK_VER}.tar.gz
tar xvzf lapack-${LAPACK_VER}.tar.gz
mkdir lib
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