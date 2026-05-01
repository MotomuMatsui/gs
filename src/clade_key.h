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

#ifndef CLADE_KEY_H
#define CLADE_KEY_H

#include <cstdint>
#include <vector>

// FNV-1a 64-bit hash of a numerically-sorted ID list. Used as clade key by
// the EP/bootstrap accumulators (EP_fbs, EP_tbe) and the lookup in addEP().
//
// Replaces the previous "1|3|7|" stringstream key:
//   * no per-clade std::string allocation
//   * no digit-by-digit base-10 conversion
//   * fixed-cost hash regardless of ID magnitude
//
// The collision probability is negligible at our scale (≤10^6 unique clades
// against 2^64 buckets ⇒ <10⁻⁷ by birthday bound).
//
// IMPORTANT: callers must pass `ids` already sorted ascending — both EP_fbs
// (iterates z=0..size-1) and addEP (calls std::sort on clade[depth]) feed in
// numerically-sorted IDs, so the requirement is met without extra work.
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

#endif /*CLADE_KEY_H*/
