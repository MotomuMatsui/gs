<p align="center"><img src="https://raw.githubusercontent.com/MotomuMatsui/gs/master/GSbanner.png"></p>  

# gs
Graph Splitting (GS) is a brand-new phylogenetic analysis method, which can effectively resolve early evolution of protein families. Its accuracy and speed was proved by extensive evolutionary simulation, and its application to TIM-barrel superfamily highlighted instant evolution of protein-mediated pyrimidine biosynthesis at a transitional phase between the RNA world and the modern DNA-RNA-protein world.

<p align="center"><img src="https://raw.githubusercontent.com/MotomuMatsui/gs/master/introduction.png"></p>

Online tool: [GS analysis server](http://gs.bs.s.u-tokyo.ac.jp/)  
Reference: Matsui and Iwasaki, ???, 2018  
Our Laboratory: [Iwasaki Lab](http://iwasakilab.bs.s.u-tokyo.ac.jp/eindex.html)  

[![Build Status](https://travis-ci.org/MotomuMatsui/gs.svg?branch=master)](https://travis-ci.org/MotomuMatsui/gs)
[![Language](https://img.shields.io/badge/C%2B%2B-5.0%2B-green.svg)](https://gcc.gnu.org/)
[![LAPACK](https://img.shields.io/badge/LAPACK%2FBLAS-3.8%2B-green.svg)](http://www.netlib.org/lapack/)
[![MMseqs](https://img.shields.io/badge/MMSeqs-2.0%2B-green.svg)](https://github.com/soedinglab/MMseqs2)
[![Ubuntu](https://img.shields.io/badge/Linux-Ubuntu-green.svg)](https://www.ubuntu.com/)
[![CentOS](https://img.shields.io/badge/Linux-CentOS-green.svg)](https://www.centos.org/)
[![Mac](https://img.shields.io/badge/Mac-macOS-green.svg)](https://www.apple.com/macos/)
[![GPL License](https://img.shields.io/badge/license-GPL-blue.svg)](LICENSE)

## Version
version 1.0 (2018/05/24)

## Installation

### Supported OS
`gs` is available for Linux and Mac (macOS).

### Requirements
1. [GNU GCC compiler](https://gcc.gnu.org/) (5.0+)
1. [MMseqs2 package](https://github.com/soedinglab/mmseqs2) (2.0+)
1. [LAPACK/BLAS package](http://www.netlib.org/lapack/) (3.8+)
  
:exclamation: LAPACK/BLAS package will be installed at the next section.

### Installation

    $ git clone https://github.com/MotomuMatsui/gs
    $ cd gs
    $ sh ./lapack_install.sh
    $ make

## Usage
    $ ./gs [-e INTEGER(>=0)] [-s (silent mode)] [-h (help)] IN(fasta) > OUT(newick)

## Example
Phylogenetic tree WITH branch reliability (Edge perturbation; EP) scores:

    $ ./gs -e 100 test/test.fst > test/test.nwk

Phylogenetic tree WITHOUT EP scores + silent mode:

    $ ./gs -e 0 -s test/test.fst > test/test.nwk

Show help messages:

    $ ./gs -h

## License
This software is distributed under the GNU GPL, see [LICENSE](LICENSE).
Copyright (c) 2018, Motomu Matsui

## Author
[Motomu Matsui](https://sites.google.com/site/motomumatsui/)

## Reference
Motomu Matsui and Wataru Iwasaki. ??? (2018)

## Acknowledgements
This package uses LAPACKE/CBLAS 3.8.0 provided by Univ. of Tennessee; Univ. of California, Berkeley; Univ. of Colorado Denver; and NAG Ltd.. The authors thank LAPACK/BLAS team. You can get the detailed information from http://www.netlib.org/lapack/.
