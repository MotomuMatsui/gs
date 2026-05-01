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

#ifndef FORMAT_H
#define FORMAT_H

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

// FNV-1a 64-bit hash of a numerically-sorted ID list.  Shared between the EP
// accumulators (EP_fbs, EP_tbe in ep.cpp) and addEP() below.  Replaces the
// former "1|3|7|" stringstream key — no allocation, no digit conversion.
// IDs MUST be passed already sorted ascending; both call sites satisfy this.
namespace clade_key {
inline std::uint64_t hash_sorted_ids(const std::vector<int>& ids){
  std::uint64_t h = 0xcbf29ce484222325ULL; // FNV-1a 64-bit offset basis
  for(int x : ids){
    h ^= static_cast<std::uint64_t>(x);
    h *= 0x100000001b3ULL;                 // FNV-1a 64-bit prime
  }
  return h;
}
} // namespace clade_key

int readFASTA(ifstream& ifs, ofstream& ofs1, ofstream& ofs2, int& row, bool& is_nucleotide);
int bl2mat(ifstream& ifs, double* (&W), int const& size);
int distmat2mat(ifstream& ifs, ofstream& ofs_ann, double* (&W), int& size);
int simmat2mat(ifstream& ifs, ofstream& ofs_ann, double* (&W), int& size);
void sc2nwk(int* const& W, string& newick, int const& size);
void addEP(string const& newick, string& newick_EP, unordered_map<std::uint64_t, double>& ep, int const& ep_num, int const& size);
void addLABEL(string const& newick, string& newick_ann, string const& annotation_txt, int const& size);
void sc2list(int* const& gs, int* (&list), int const& size);

#endif /*FORMAT_H*/
