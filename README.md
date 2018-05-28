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

## Requirements
- `gs1` is available for Linux and Mac (macOS).

- You will need the following tools to complie `gs1`:
    1. [GNU GCC compiler](https://gcc.gnu.org/) (5.0+)
    1. [MMseqs2](https://github.com/soedinglab/mmseqs2) (2.0+)
    1. [LAPACK/BLAS package](http://www.netlib.org/lapack/) (3.8+)

:exclamation: If you want to compile `gs1` on Mac, please install `gcc` from [Homebrew](https://brew.sh/).  
:exclamation: LAPACK/BLAS package will be installed in the next section.

## Installation

    $ git clone https://github.com/MotomuMatsui/gs
    $ cd gs
    $ sh ./lapack_install.sh
    $ make

## Usage
To get on-line help:

    $ ./gs1 -h
    
The following command enables you to calculate GS tree (phylogenetic tree reconstructed by Graph Splitting method):

    $ ./gs1 [-e INTEGER(>=0)] [-r INTEGER(>0)] [-t INTEGAR(>0)] [-s] [-h] [-v] input > output

:exclamation: A multiple sequence file (ex. test/test.fst) is required in [fasta format](https://en.wikipedia.org/wiki/FASTA_format).

## Example
GS tree (in [newick format](https://en.wikipedia.org/wiki/Newick_format)) will be displayed in `STDOUT`:

    $ ./gs1 test/test.fst

GS tree with branch reliability (Edge perturbation; EP) scores will be saved in `test.nwk`:

    $ ./gs1 -e 100 test/test.fst > test/test.nwk

GS tree with EP scores, a seed number is specified for EP method:

    $ ./gs1 -e 100 -r 12345 test/test.fst > test/test.nwk

GS tree WITHOUT EP scores + silent mode:

    $ ./gs1 -e 0 -s test/test.fst > test/test.nwk

MMseqs2 runs multithreaded jobs (4 CPUs are used in parallel):

    $ ./gs1 -e 100 -t 4 test/test.fst > test/test.nwk

Show help messages:

    $ ./gs1 -h

Show the version of `gs1`:

    $ ./gs1 -v

## License
This software is distributed under the GNU GPL, see [LICENSE](LICENSE).  
Copyright &copy; 2018, Motomu Matsui

## Author
[Motomu Matsui](https://sites.google.com/site/motomumatsui/)

## Reference
Motomu Matsui and Wataru Iwasaki. ??? (2018)

## Acknowledgements
This package includes the LAPACKE/CBLAS package provided by Univ. of Tennessee; Univ. of California, Berkeley; Univ. of Colorado Denver; and NAG Ltd.. The authors thank the LAPACK/BLAS team. You can get the detailed information from http://www.netlib.org/lapack/.
