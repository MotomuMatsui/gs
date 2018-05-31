<p align="center"><img src="https://raw.githubusercontent.com/MotomuMatsui/gs/master/GSbanner.png"></p>  

# gs2
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


## Demo

![demo](https://raw.githubusercontent.com/MotomuMatsui/gs/master/demo.gif)

## Version
version 2.0 (2018/05/30)

## Requirements

- `gs2` is available for Linux and Mac (macOS).
- [GNU GCC compiler](https://gcc.gnu.org/) (5.0+) is required to compile `gs2`

:exclamation: If you want to compile `gs2` on Mac, please install `gcc` from [Homebrew](https://brew.sh/).  

## Installation on Linux/Mac from source code

#### 1. Obtain the gs package:

    $ git clone https://github.com/MotomuMatsui/gs
    $ cd gs

#### 2. Install the MMseqs2:

If [MMseqs2](https://github.com/soedinglab/mmseqs2) (2.0+) has been already installed in your system, you can skip this step.

    $ sh ./mmseqs_install.sh
    $ export PATH=$(pwd)/MMseqs2/build/bin:$PATH

:exclamation: If this script made some errors, please install `mmseqs` with reference to [MMseqs page](https://github.com/soedinglab/mmseqs2).   
:exclamation: Optionally, you can move `mmseqs` to the other place where you want (ex. `~/bin`) and add this path to your path environment variable (ex. `export PATH=~/bin:$PATH`).

#### 3. Install the LAPACK/BLAS package:

If [LAPACK/BLAS package](http://www.netlib.org/lapack/) (3.8+) package has been already installed in your system, you can skip this step.

    $ sh ./lapack_install.sh

:exclamation: If this script made some errors, please check the [LAPACK/BLAS page](http://www.netlib.org/lapack).

#### 4. Install the gs package:

    $ make

:exclamation: If necessary, please modify the Makefile in response to your environment (ex. `CXX := g++-8`, `INC := -I/usr/local/include`, `CXXFLAGS += -std=c++1z`).

### Known issues

- If you have some errors when compiling LAPACK/BLAS pakage on MacOS, please rewrite `OPTS = -O2 -frecursive` to `OPTS = -O3 -frecursive -pipe` in `lapack-3.8.0/make.inc`, then re-execute `source ./lapack_install.sh` and `make`.

## Usage
To get on-line help:

    $ ./gs2 -h
    
The following command enables you to calculate GS tree (phylogenetic tree reconstructed by Graph Splitting method):

    $ ./gs2 [arguments] input > output

:exclamation: A multiple sequence file (ex. [data/200.faa](data/200.faa)) should be required as `input` in [fasta format](https://en.wikipedia.org/wiki/FASTA_format).

### Arguments

|Option| Description                                                                                         |
|:----:|:----------------------------------------------------------------------------------------------------|
|  -e  |<strong>[integer(>=0)]</strong> <em>The number of replicates for EP method. Default: 0</em>          |
|  -r  |<strong>[integer(>=1)]</strong> <em>The random seed number for EP method. Default: random number</em>|
|  -t  |<strong>[integer(>=1)]</strong> <em>The number of threads for MMseqs. Default: 1</em>                |
|  -m  |<strong>[float(1&ndash;7.5)]</strong> <em>Sensitivity for MMseqs. Default: 7.5</em>                  |
|  -s  |<em>Silent mode: do not report progress. Default: Off</em>                                           |
|  -h  |<em>Show help messages. Default: Off</em>                                                            |
|  -v  |<em>Show the version. Default: Off</em>                                                              |

## Examples
GS tree (in [newick format](https://en.wikipedia.org/wiki/Newick_format)) will be displayed in `STDOUT`:

    $ ./gs2 data/200.faa

GS tree with branch reliability (Edge perturbation; EP) scores will be saved in `test.nwk`:

    $ ./gs2 -e 100 data/200.faa > data/200.nwk

GS tree with EP scores; a seed number is specified for EP method:

    $ ./gs2 -e 100 -r 12345 data/200.faa > data/200.nwk

GS tree WITHOUT EP scores + silent mode:

    $ ./gs2 -e 0 -s data/200.faa > data/200.nwk

MMseqs2 runs multithreaded jobs (4 CPUs are used in parallel):

    $ ./gs2 -e 100 -t 4 data/200.faa > data/200.nwk

Visualization of [200.nwk](data/200.nwk) by [iTOL](https://itol.embl.de/):

<p align="center"><img src="https://raw.githubusercontent.com/MotomuMatsui/gs/master/data/200_iTOL.png"></p>  

## License
This software is distributed under the GNU GPL, see [LICENSE](LICENSE).  
Copyright &copy; 2018, Motomu Matsui

## Author
[Motomu Matsui](https://sites.google.com/site/motomumatsui/)

## Reference
Motomu Matsui and Wataru Iwasaki. ??? (2018)

## Acknowledgements
This package includes the LAPACKE/CBLAS (Univ. of Tennessee; Univ. of California, Berkeley; Univ. of Colorado Denver; and NAG Ltd.) and MMseqs (S&ouml;ding Laboratory) packages. The authors give special thanks to both teams. You can get the detailed information from http://www.netlib.org/lapack/ and https://github.com/soedinglab/MMseqs2.
