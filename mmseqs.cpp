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

#include <stdio.h>
#include <fstream>
#include <sstream>
#include <string>
#include <regex>
#include <unordered_map>
#include <sys/stat.h>

using namespace std;

void mmseqs(string const& fst, string const& output){
  //mmseqs command
  const string mmseqs = "mmseqs";

  // File I/O
  regex re(R"(_simple\.fst$)"); // Matching to the extention of the input file
  auto queryDB  = regex_replace(fst, re, "_queryDB");
  auto targetDB = regex_replace(fst, re, "_targetDB");
  auto resultDB = regex_replace(fst, re, "_resultDB");
  auto tmp      = regex_replace(fst, re, "_tmp");

  //Commands
  auto cmd1 = mmseqs+" createdb "    +fst+" "+queryDB+" > /dev/null";
  auto cmd2 = mmseqs+" createdb "    +fst+" "+targetDB+" > /dev/null";
  auto cmd3 = mmseqs+" search "      +queryDB+" "+targetDB+" "+resultDB+" "+tmp+" --threads 1 -e 10 -s 7.5 > /dev/null";
  auto cmd4 = mmseqs+" convertalis " +queryDB+" "+targetDB+" "+resultDB+" "+output+" --format-mode 0 > /dev/null";

  mkdir(tmp.c_str(), S_IRWXU | S_IRGRP | S_IROTH);
  system(cmd1.c_str()); //MMseqs: createdb
  system(cmd2.c_str()); //MMseqs: createdb
  system(cmd3.c_str()); //MMseqs: search
  system(cmd4.c_str()); //MMseqs: convertalis

  //Delete intermediate files
  string extension[6] = {
    "", ".dbtype", ".lookup",
    ".index", "_h.index", "_h"};
  
  for(auto e : extension){
    auto file1 = queryDB  + e;
    auto file2 = targetDB + e;
    auto file3 = resultDB + e;
    remove(file1.c_str());
    remove(file2.c_str());
    remove(file3.c_str());
  }
}
