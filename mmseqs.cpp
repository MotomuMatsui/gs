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

#include <stdio.h>
#include <fstream>
#include <sstream>
#include <string>
#include <regex>
#include <unordered_map>
#include <sys/stat.h>

using namespace std;

void mmseqs(string const& fst, string const& output, string const& threads, string const& sensitivity){
  //mmseqs command
  const string mmseqs = "mmseqs";

  // File I/O
  regex re(R"(_simple\.fst$)"); // Matching to the extention of the input file
  auto DB     = regex_replace(fst, re, "_DB");
  auto result = regex_replace(fst, re, "_result");
  auto tmp    = regex_replace(fst, re, "_tmp");

  //Commands
  auto cmd1 = mmseqs+" createdb "    +fst+" "+DB+" > /dev/null";
  auto cmd2 = mmseqs+" search "      +DB+" "+DB+" "+result+" "+tmp+" --threads "+threads+" -e 10 -s "+sensitivity+" > /dev/null";
  auto cmd3 = mmseqs+" convertalis " +DB+" "+DB+" "+result+" "+output+" --format-mode 0 > /dev/null";

  mkdir(tmp.c_str(), S_IRWXU | S_IRGRP | S_IROTH); //mkdir tmp
  auto info1 = system(cmd1.c_str()); if(info1>0){} //MMseqs: createdb
  auto info2 = system(cmd2.c_str()); if(info2>0){} //MMseqs: search
  auto info3 = system(cmd3.c_str()); if(info3>0){} //MMseqs: convertalis

  //Delete intermediate files
  string extension[6] = {
    "", ".dbtype", ".lookup",
    ".index", "_h.index", "_h"};
  
  for(auto e : extension){
    auto file1 = DB     + e;
    auto file2 = result + e;
    remove(file1.c_str());
    remove(file2.c_str());
  }
}
