/*=========================================*\
| Graph Splitting Method v1.0 (2018/06/11) |
|                                          |
| Copyright (c) 2018 Motomu Matsui         |
| Systematic Biology, xx:xx-xx, 2018       |
|                                          |
| Web:  http://gs.bs.s.u-tokyo.ac.jp/      |
| Mail: qm at bs.s.u-tokyo.ac.jp           |
\*=========================================*/

#include <iostream>
#include <string>

using namespace std;

void print_banner(){

  string banner = 
    "--------------------------------------------------\n"
    " Graph Splitting Method v1.0 (2018/06/01)         \n"
    " http://gs.bs.s.u-tokyo.ac.jp/                    \n"
    " Matsui M and Iwasaki W, Systematic Biology, 2018 \n"
    "--------------------------------------------------\n";

  cerr << banner << endl;
}

void print_usage(char*& program){
  cerr << "Usage: " << program << " [-e INTEGER(>=0)] [-s (silent mode)] [-h (help)] IN(fasta) > OUT(newick)" << endl;
}
