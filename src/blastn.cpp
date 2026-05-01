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

#include <fstream>
#include <regex>
#include <sstream>
#include <stdio.h>
#include <string>
#include <sys/stat.h>
#include <unordered_map>

#include "blastn.h"

using namespace std;

void blastn(string const& fst, string const& output, string const& threads){
  // BLAST+ commands
  const string makeblastdb = "makeblastdb";
  const string blastn_cmd  = "blastn";

  // File I/O
  regex re(R"(_simple\.fst$)"); // Matching the extension of the input file
  auto DB = regex_replace(fst, re, "_blastdb");

  // Commands
  // makeblastdb builds a nucleotide database from the simplified FASTA.
  // blastn runs an all-vs-all search; -outfmt 6 produces the same 12-column
  // tabular layout that bl2mat() already parses (query, target, ..., bitscore).
  auto cmd1 = makeblastdb + " -in " + fst + " -dbtype nucl -out " + DB + " > /dev/null 2>&1";
  auto cmd2 = blastn_cmd  + " -query " + fst + " -db " + DB
              + " -outfmt 6 -evalue 10"
              + " -num_threads " + threads
              + " -out " + output + " > /dev/null 2>&1";

  auto info1 = system(cmd1.c_str()); if(info1>0){} // makeblastdb
  auto info2 = system(cmd2.c_str()); if(info2>0){} // blastn

  // Delete intermediate database files produced by makeblastdb.
  string extension[10] = {
    ".nhr", ".nin", ".nsq",
    ".ndb", ".nhd", ".nhi",
    ".not", ".nto", ".ntf",
    ".njs"};

  for(auto e : extension){
    auto file = DB + e;
    remove(file.c_str());
  }
}
