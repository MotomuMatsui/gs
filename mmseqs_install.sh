#!/usr/bin/sh

# Decompress MMseqs package
tar xzf MMseqs2.tar.gz
cd MMseqs2
if [ ! -d build ]; then mkdir build; fi
cd build

# CMAKE
cmake -DHAVE_SSE4_1=1 -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..

# Compile MMseqs
make
make install 

cd ..
cd ..
