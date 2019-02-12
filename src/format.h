/******************************************\
| Graph Splitting Method v2.4 (2019/02/12) |
|                                          |
| Copyright (c) 2015-2019 Motomu Matsui    |
|     Distributed under the GNU GPL        |
|                                          |
|     Matsui M and Iwasaki W (2019)        |
|     Systematic Biology, xx:xx-xx.        |
|                                          |
|     http://gs.bs.s.u-tokyo.ac.jp/        |
\******************************************/

#ifndef FORAMT_H
#define FORMAT_H

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

int readFASTA(ifstream& ifs, ofstream& ofs1, ofstream& ofs2, int& row);
int bl2mat(ifstream& ifs, double* (&W), int const& size);
void sc2nwk(int* const& W, string& newick, int const& size);
void addEP(string const& newick, string& newick_EP, unordered_map<string, double>& ep, int const& ep_num, int const& size);
void addLABEL(string const& newick, string& newick_ann, string const& annotation_txt, int const& size);
void sc2list(int* const& gs, int* (&list), int const& size);

#endif /*FORMAT_H*/
