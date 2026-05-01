/******************************************\
| Graph Splitting Method v2.5 (2026/04/23) |
|                                          |
| Copyright (c) Motomu Matsui              |
|     Distributed under the GNU GPL        |
|                                          |
|     Matsui M and Iwasaki W (2019)        |
|     Systematic Biology, 69(2):265-279.   |
|                                          |
|     http://gs.bs.s.u-tokyo.ac.jp/        |
\******************************************/

#ifndef BLASTN_H
#define BLASTN_H

#include <fstream>
#include <regex>
#include <sstream>
#include <stdio.h>
#include <string>
#include <sys/stat.h>
#include <unordered_map>

using namespace std;

void blastn(string const& fst, string const& output, string const& threads);

#endif /*BLASTN_H*/
