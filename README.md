# gs
Graph Splitting (GS) is a brand-new phylogenetic analysis method, which can effectively resolve early evolution of protein families. Its accuracy and speed was proved by extensive evolutionary simulation, and its application to TIM-barrel superfamily highlighted instant evolution of protein-mediated pyrimidine biosynthesis at a transitional phase between the RNA world and the modern DNA-RNA-protein world.

Online tool: [GS analysis server](http://gs.bs.s.u-tokyo.ac.jp/)

## Version
version 1.0 (2018/05/24)

## Acknowledgements
The gs package uses LAPACKE/CBLAS 3.8.0 provided by Univ. of Tennessee; Univ. of California, Berkeley; Univ. of Colorado Denver; and NAG Ltd.. The authors thank LAPACK/BLAS team. You can get the detailed information from http://www.netlib.org/lapack/ and http://www.netlib.org/blas/.

## Installation

    git clone https://github.com/MotomuMatsui/gs
    cd gs
    mkdir lib
    sh ./lapack_install.sh
    make
    cp gs /your/favorite/path/

## Usage
    ./gs [-e INTEGER(>=0)] [-s (silent mode)] [-h (help)] IN(fasta) > OUT(newick)

## Example
    Phylogenetic tree with branch reliability scores (EP values):
    ./gs -e 100 test.fst > test.nwk  
    
    Only topology of phylogenetic tree + silent mode:
    ./gs -e 0 -s test.fst > test.nwk
    
    Show help messages:
    ./gs -h

## License
This software is distributed under the GNU GPL, see LICENSE.

## Reference
Motomu Matsui and Wataru Iwasaki. ??? (2018)
