/******************************************\
| Graph Splitting Method v1.0 (2018/06/11) |
|                                          |
| Copyright (c) 2015-2018 Motomu Matsui    |
|     Distributed under the GNU GPL        |
|                                          |
|     Matsui M and Iwasaki W (2018)        |
|     Systematic Biology, xx:xx-xx.        |
|                                          |
|     http://gs.bs.s.u-tokyo.ac.jp/        |
\******************************************/

#include <iostream>
#include <string>

using namespace std;

void print_banner(){

  string banner = 
    "------------------------------------------\n"
    " Graph Splitting Method v1.0 (2018/06/01) \n"
    "                                          \n"
    "   Copyright (c) 2018 Motomu Matsui       \n"
    "   Systematic Biology, xx:xx-xxx, 2018    \n"
    "                                          \n"
    "   http://gs.bs.s.u-tokyo.ac.jp/          \n"
    "------------------------------------------\n";

  cerr << banner << endl;
}

void print_usage(char*& program){
  cerr << "Usage: " << program << " [-e INTEGER(>=0)] [-s (silent mode)] [-h (help)] IN(fasta) > OUT(newick)" << endl;
}
