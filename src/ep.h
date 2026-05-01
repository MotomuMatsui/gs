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

#ifndef EP_H
#define EP_H

#include <cstdint>
#include <unordered_map>

// seed:   RNG seed (each OpenMP thread is offset by thread_id * 1000003)
// ep_num: number of EP iterations (parallelised internally)
// ep:     bootstrap accumulator keyed by clade hash (see clade_key.h)
void EP_fbs(double* const (&W), std::unordered_map<std::uint64_t,double>& ep,
            int const seed, int const ep_num, int const& size);

void EP_tbe(double* const (&W), int* const (&list_ori),
            std::unordered_map<std::uint64_t,double>& ep,
            int const seed, int const ep_num, int const& size);

#endif /*EP_H*/
