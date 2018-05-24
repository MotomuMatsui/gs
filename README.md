# gs
Graph Splitting (GS) is a brand-new phylogenetic analysis method, which can effectively resolve early evolution of protein families. Its accuracy and speed was proved by extensive evolutionary simulation, and its application to TIM-barrel superfamily highlighted instant evolution of protein-mediated pyrimidine biosynthesis at a transitional phase between the RNA world and the modern DNA-RNA-protein world.

## Version
version 1.0 (2018/05/24)

## Acknowledgements
The gs package includes LAPACK and CBLAS 3.8.0 provided by Univ. of Tennessee; Univ. of California, Berkeley; Univ. of Colorado Denver; and NAG Ltd..
You can get the detailed information from http://www.netlib.org/lapack/ and http://www.netlib.org/blas/.
The authors thank LAPACK/BLAS team.

## Usage
./gs [-e INTEGER(>=0)] [-s (silent mode)] [-h (help)] IN(fasta) > OUT(newick)

## Example
./gs -e 100 test.fst > test.nwk  # Phylogenetic tree with branch reliability scores (EP values)
./gs -e 0 -s test.fst > test.nwk # Only topology of phylogenetic tree + silent mode
./gs -h                          # Show help messages

## License
This software is distributed under the GNU GPL, see LICENSE.

## Reference
Motomu Matsui and Wataru Iwasaki. ??? (2018)
