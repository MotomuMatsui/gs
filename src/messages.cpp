/******************************************\
| Graph Splitting Method v2.0 (2018/06/01) |
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
    " Graph Splitting Method v2.0 (2018/06/01) \n"
    "                                          \n"
    "   Copyright (c) 2018 Motomu Matsui       \n"
    "   Systematic Biology, xx:xx-xxx, 2018    \n"
    "                                          \n"
    "   http://gs.bs.s.u-tokyo.ac.jp/          \n"
    "------------------------------------------\n";

  cerr << banner << endl;
}

void print_usage(char*& program){
  cerr << "Usage: " << program << " [-e INTEGER(>=0)] [-r INTEGER(>0)] [-t INTEGAR(>0)] [-s] [-h] [-v] input > output" << endl;
  cerr << "-e " << "the number of replicates for EP method. Default: 0" << endl;
  cerr << "-r " << "the random seed number for EP method. Default: random number" << endl;
  cerr << "-t " << "the number of threads for MMseqs. Default: 1" << endl;
  cerr << "-m " << "sensitivity for MMseqs. Default: 7.5" << endl;
  cerr << "-s " << "silent mode: do not report progress. Default: Off" << endl;
  cerr << "-h " << "show help messages. Default: Off" << endl;
  cerr << "-v " << "show the version. Default: Off" << endl;
}
