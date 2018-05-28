# Decompression LAPACK package.
git clone https://github.com/soedinglab/MMseqs2.git
cd MMseqs2
mkdir build
cd build

# CMAKE
cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..

# Compile MMseqs
make
make install 

# Environment variable
export PATH=$(pwd)/bin/:$PATH
cd ..
cd ..
