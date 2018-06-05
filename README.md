<p align="center"><img src="https://raw.github.com/wiki/MotomuMatsui/gs/images/GSbanner.png"></p>  

# gs2
`gs2` is a software to conduct a brand-new phylogenetic analysis method--the Graph Splitting (GS). It can effectively resolve early evolution of protein families and its accuracy and speed was proved by extensive evolutionary simulation.    
`gs2` is open-source software (GPL v3.0) inplemented in C++ for Linux, Mac and Windows.

<p align="center"><img src="https://raw.github.com/wiki/MotomuMatsui/gs/images/introduction.png"></p>

Reference: Matsui and Iwasaki, ???, 2018  
Online tool: [GS analysis server](http://gs.bs.s.u-tokyo.ac.jp/)  
Our Laboratory: [Iwasaki Lab](http://iwasakilab.bs.s.u-tokyo.ac.jp/eindex.html)  

[![Build Status](https://travis-ci.org/MotomuMatsui/gs.svg?branch=master)](https://travis-ci.org/MotomuMatsui/gs)
[![Ubuntu](https://img.shields.io/badge/Linux-Ubuntu-green.svg)](https://www.ubuntu.com/)
[![CentOS](https://img.shields.io/badge/Linux-CentOS-green.svg)](https://www.centos.org/)
[![Mac](https://img.shields.io/badge/Mac-macOS-green.svg)](https://www.apple.com/macos/)
[![Windows](https://img.shields.io/badge/Windows-Cygwin-green.svg)](https://www.cygwin.com/)
[![Language](https://img.shields.io/badge/C%2B%2B-5.0%2B-green.svg)](https://gcc.gnu.org/)
[![LAPACK](https://img.shields.io/badge/LAPACK%2FBLAS-3.7%2B-green.svg)](http://www.netlib.org/lapack/)
[![MMseqs](https://img.shields.io/badge/MMSeqs-2.0%2B-green.svg)](https://github.com/soedinglab/MMseqs2)
[![GPL License](https://img.shields.io/badge/license-GPL-blue.svg)](LICENSE)

## History
version 2.0 (2018/06/01)   
  - Re-implemented in C++    
  - MMseq2 is used for all-to-all pairwise sequence alignment    

version 1.0 (2017/02/07)   
  - Implemented in R and Perl    
  - BLAST+ is used for all-to-all pairwise sequence alignment    

## Demo

![demo](https://raw.github.com/wiki/MotomuMatsui/gs/images/demo.gif)

## Installation

### 0. Requirements

- `gs2` is available for <strong>Linux</strong>, <strong>Mac (macOS)</strong> and <strong>Windows (Cygwin)</strong>.
- [GNU GCC compiler](https://gcc.gnu.org/) (5.0+) is required to compile `gs2`

:exclamation: For Mac users: please install `gcc` from [Homebrew](https://brew.sh/).  

### 1. Download gs package:

    $ git clone https://github.com/MotomuMatsui/gs
    $ cd gs

### 2. Install MMseqs2:

    $ sh ./mmseqs_install.sh
    $ export PATH=$(pwd)/MMseqs2/build/bin:$PATH

- Optionally, you can move `mmseqs` to the other place where you want (ex. `~/bin`) and add this path to your path environment variable (ex. `export PATH=~/bin:$PATH`).

### 3. Install LAPACK/BLAS package:

    $ sh ./lapack_install.sh

### 4. Install gs package:

    $ make

- If necessary, you can optimize the Makefile in response to your environment (ex. `CXX := g++-8`, `CXXFLAGS += -std=c++1z`).

### Known issues

- For Mac users: if you have some errors when compiling LAPACK/BLAS pakage, please rewrite `OPTS = -O2 -frecursive` to `OPTS = -O3 -frecursive -pipe` in `lapack-3.7.1/make.inc`, then re-execute `source ./lapack_install.sh` and `make`.    
- For Windows users: LAPACK/BLAS <strong>version 3.8.0</strong> has some problem to be installed into Windows. Please choose LAPACK/BLAS <strong>version 3.7.1</strong> (default) for installation.

## Usage
To get on-line help:

    $ ./gs2 -h
    
The following command enables you to calculate GS tree (phylogenetic tree reconstructed by Graph Splitting method):

    $ ./gs2 [arguments] input > output

:exclamation: A multiple sequence file (ex. [example/200.faa](example/200.faa)) should be required as `input` in [fasta format](https://en.wikipedia.org/wiki/FASTA_format).

Arguments:

|Option| Description                                                                                         |
|:----:|:----------------------------------------------------------------------------------------------------|
|  -e  |<strong>[integer(>=0)]</strong> <em>The number of replicates for EP method. Default: 0</em>          |
|  -r  |<strong>[integer(>=1)]</strong> <em>The random seed number for EP method. Default: random number</em>|
|  -t  |<strong>[integer(>=1)]</strong> <em>The number of threads for MMseqs. Default: 1</em>                |
|  -m  |<strong>[float(1&ndash;7.5)]</strong> <em>Sensitivity for MMseqs. Default: 7.5</em>                  |
|  -s  |<em>Silent mode: do not report progress. Default: Off</em>                                           |
|  -l  |<em>use actual names in newick format. Default: Off</em>                                             |
|  -h  |<em>Show help messages. Default: Off</em>                                                            |
|  -v  |<em>Show the version. Default: Off</em>                                                              |

## Examples
GS tree (in [newick format](https://en.wikipedia.org/wiki/Newick_format)) will be displayed in `STDOUT`:

    $ ./gs2 example/200.faa

GS tree with branch reliability (Edge perturbation; EP) scores will be saved in `test.nwk`:

    $ ./gs2 -e 100 example/200.faa > example/200.nwk

GS tree with EP scores; a seed number is specified for EP method:

    $ ./gs2 -e 100 -r 12345 example/200.faa > example/200.nwk

GS tree WITHOUT EP scores + silent mode:

    $ ./gs2 -e 0 -s example/200.faa > example/200.nwk

MMseqs2 runs multithreaded jobs (4 CPUs are used in parallel):

    $ ./gs2 -e 100 -t 4 example/200.faa > example/200.nwk

Visualization of [200.nwk](example/200.nwk) by [iTOL](https://itol.embl.de/):

<p align="center"><img src="https://raw.github.com/wiki/MotomuMatsui/gs/images/200_iTOL.png"></p>  

## License
This software is distributed under the GNU GPL, see [LICENSE](LICENSE).  
Copyright &copy; 2018, Motomu Matsui

## Author
[Motomu Matsui](https://sites.google.com/site/motomumatsui/)

## Reference
Motomu Matsui and Wataru Iwasaki. ??? (2018)

## Acknowledgements
This package includes the LAPACKE/CBLAS (Univ. of Tennessee; Univ. of California, Berkeley; Univ. of Colorado Denver; and NAG Ltd.) and MMseqs (S&ouml;ding Laboratory) packages. The authors give special thanks to both teams. You can get the detailed information from http://www.netlib.org/lapack/ and https://github.com/soedinglab/MMseqs2.
