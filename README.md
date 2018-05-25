*****
<p align="center"><img src="https://raw.githubusercontent.com/MotomuMatsui/gs/master/GSbanner.png"></p>
*****

# gs
Graph Splitting (GS) is a brand-new phylogenetic analysis method, which can effectively resolve early evolution of protein families. Its accuracy and speed was proved by extensive evolutionary simulation, and its application to TIM-barrel superfamily highlighted instant evolution of protein-mediated pyrimidine biosynthesis at a transitional phase between the RNA world and the modern DNA-RNA-protein world.

Online tool: [GS analysis server](http://gs.bs.s.u-tokyo.ac.jp/)

[![Build Status](https://travis-ci.org/MotomuMatsui/gs.svg?branch=master)](https://travis-ci.org/MotomuMatsui/gs)
[![Language](https://img.shields.io/badge/C%2B%2B-%E2%89%A55.0-green.svg)](https://gcc.gnu.org/)
[![LAPACK](https://img.shields.io/badge/LAPACK%2FBLAS-%E2%89%A53.8-green.svg)](http://www.netlib.org/lapack/)
[![MMseqs](https://img.shields.io/badge/MMSeqs-%E2%89%A52-green.svg)](https://github.com/soedinglab/MMseqs2)
[![Ubuntu](https://img.shields.io/badge/Linux-Ubuntu-green.svg)](https://www.ubuntu.com/)
[![CentOS](https://img.shields.io/badge/Linux-CentOS-green.svg)](https://www.centos.org/)
[![Mac](https://img.shields.io/badge/Mac-macOS-green.svg)](https://www.apple.com/macos/)
[![GPL License](https://img.shields.io/badge/license-GPL-blue.svg)](LICENSE)

## Version
version 1.0 (2018/05/24)

## Acknowledgements
This package uses LAPACKE/CBLAS 3.8.0 provided by Univ. of Tennessee; Univ. of California, Berkeley; Univ. of Colorado Denver; and NAG Ltd.. The authors thank LAPACK/BLAS team. You can get the detailed information from http://www.netlib.org/lapack/.

## Installation

    git clone https://github.com/MotomuMatsui/gs
    cd gs
    mkdir lib
    sh ./lapack_install.sh
    make
    cp gs /your/favorite/path/
    export PATH=/your/favorite/path/:$PATH

## Usage
    gs [-e INTEGER(>=0)] [-s (silent mode)] [-h (help)] IN(fasta) > OUT(newick)

## Example
Phylogenetic tree WITH branch reliability (Edge perturbation; EP) scores:

    gs -e 100 test.fst > test.nwk

Phylogenetic tree WITHOUT EP scores + silent mode:

    gs -e 0 -s test.fst > test.nwk

Show help messages:

    gs -h

## License
This software is distributed under the GNU GPL, see LICENSE.

## Reference
Motomu Matsui and Wataru Iwasaki. ??? (2018)
