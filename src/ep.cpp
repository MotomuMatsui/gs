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

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iostream>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#ifdef _OPENMP
#  include <omp.h>
#endif

#include "sc.h"
#include "gs_functions.h"
#include "format.h" // pulls in clade_key::hash_sorted_ids

using namespace std;

// Run one GS iteration on perturbed matrix E; returns list_ep ownership.
// Shared by EP_fbs and EP_tbe so the inner loop is not duplicated.
static unique_ptr<int[]> run_gs_list(double* E, int const size){

  int N = size*size;
  vector<double> simL(size, 3);
  vector<int>    res(size, 1);

  auto step_buf = unique_ptr<int[]>(new int[N]());
  int* step_raw = step_buf.get();
  sedMATRIX(step_raw, res, res, 1, size);

  int gK   = 1;
  int gMin = 1;
  vector<int> a, b;

  while(gK < size){
    tie(a,b) = spectral_clustering(E, res, gMin);
    gK++;
    sedMATRIX(step_raw, res, b, gK, size);
    simL[gMin-1] = simI(E, res, gMin);
    simL[gK-1]   = simI(E, res, gK);
    gMin = whichMIN(simL);
  }

  int* list_ep_raw = nullptr;
  sc2list(step_raw, list_ep_raw, size);
  return unique_ptr<int[]>(list_ep_raw); // transfer ownership to caller
}

//////////////////////////////////////////////////////////////////
// Perturb similarity scores for one EP iteration.
static unique_ptr<double[]> perturb(double* const W, int const N,
                                    uniform_real_distribution<double>& urd,
                                    mt19937_64& mt){
  // Every E[i] is assigned in the loop below, so no need to zero-init.
  auto E_buf = unique_ptr<double[]>(new double[N]);
  double* E  = E_buf.get();
  for(int i = 0; i < N; i++){
    double a = W[i];
    if(a == 0){
      E[i] = 0; // Graph topology is not changed
    }
    else{
      double b = gev(urd(mt), a);
      E[i] = (b>1)? 1: (b<0)? 0: b;
    }
  }
  return E_buf;
}

//////////////////////////////////////////////////////////////////
// EP_fbs: Felsenstein bootstrap support via Edge Perturbation
// Parallelised: each OpenMP thread accumulates into a local map,
// which is merged into ep under a critical section at the end.
//////////////////////////////////////////////////////////////////

void EP_fbs(double* const (&W), unordered_map<std::uint64_t,double>& ep,
            int const seed, int const ep_num, int const& size){

  int N = size*size;

  #pragma omp parallel default(none) shared(W, ep, seed, ep_num, size, N)
  {
    unordered_map<std::uint64_t,double> local_ep;

    // Per-thread RNG — prime offset avoids correlated sequences across threads
    #ifdef _OPENMP
    mt19937_64 mt(static_cast<uint64_t>(seed) + omp_get_thread_num() * 1000003ULL);
    #else
    mt19937_64 mt(static_cast<uint64_t>(seed));
    #endif
    uniform_real_distribution<double> urd(0, 1);

    // Reusable per-thread ID buffers — avoid allocating two vectors per clade.
    std::vector<int> ids_in, ids_out;
    ids_in.reserve(size);
    ids_out.reserve(size);

    // chunk=4: each iteration runs a full GS, so per-iteration work is very
    // large and roughly uniform — the default chunk=1 wastes scheduler
    // overhead on each handout.  4 keeps load-balance flexibility while
    // amortising the OpenMP runtime cost.
    #pragma omp for schedule(dynamic, 4)
    for(int n = 0; n < ep_num; n++){
      auto E_buf   = perturb(W, N, urd, mt);
      auto list_ep = run_gs_list(E_buf.get(), size);

      for(int x = 0; x < size-3; x++){
        ids_in.clear();
        ids_out.clear();
        // z iterates 0..size-1, so both lists are filled in ascending order
        // — exactly the input precondition of clade_key::hash_sorted_ids().
        int* list_row = list_ep.get() + static_cast<size_t>(x) * size;
        for(int z = 0; z < size; z++){
          if(list_row[z] > 0) ids_in.push_back(z+1);
          else                ids_out.push_back(z+1);
        }
        local_ep[clade_key::hash_sorted_ids(ids_in)]++;
        local_ep[clade_key::hash_sorted_ids(ids_out)]++;
      }
    }

    #pragma omp critical
    {
      for(auto const& kv : local_ep) ep[kv.first] += kv.second;
    }
  } // end parallel
}

//////////////////////////////////////////////////////////////////
// EP_tbe: Transfer Bootstrap Expectation via Edge Perturbation
//////////////////////////////////////////////////////////////////

// Pack a row-major 0/1 matrix of shape (rows × size) into a bitset
// of `rows` rows × `words` 64-bit words. Returns popcount of each row.
// `bits` and `pop` must be sized rows*words and rows respectively.
static void pack_bits(int* const src, int const rows, int const size,
                      int const words,
                      uint64_t* bits, int* pop){
  std::fill(bits, bits + static_cast<size_t>(rows) * words, 0ULL);
  for(int r = 0; r < rows; r++){
    int pc = 0;
    int* src_row     = src  + static_cast<size_t>(r) * size;
    uint64_t* bit_row = bits + static_cast<size_t>(r) * words;
    for(int z = 0; z < size; z++){
      if(src_row[z] > 0){
        bit_row[z >> 6] |= 1ULL << (z & 63);
        pc++;
      }
    }
    pop[r] = pc;
  }
}

void EP_tbe(double* const (&W), int* const (&list_ori),
            unordered_map<std::uint64_t,double>& ep,
            int const seed, int const ep_num, int const& size){

  int N = size*size;

  // Bitset layout: each branch row of length `size` is packed into `words`
  // 64-bit words. Hamming distance becomes Σ popcount(ori_w ^ ep_w) — 64×
  // less memory traffic and one HW popcnt per word vs. `size` byte-level XORs
  // per pair.  Tail bits beyond `size` stay zero in both bitsets, so XOR
  // never produces spurious 1s in the unused tail.
  int const rows  = size - 3;
  int const words = (size + 63) / 64;

  // Pack list_ori once outside the parallel region (it is shared & immutable).
  std::vector<uint64_t> ori_bits(static_cast<size_t>(rows) * words, 0ULL);
  std::vector<int>      ori_pop(rows, 0);
  if(rows > 0){
    pack_bits(list_ori, rows, size, words, ori_bits.data(), ori_pop.data());
  }

  #pragma omp parallel default(none) \
          shared(W, list_ori, ep, seed, ep_num, size, N, rows, words, ori_bits, ori_pop)
  {
    unordered_map<std::uint64_t,double> local_ep;

    // Per-thread reusable bitset buffer for list_ep — avoids reallocating
    // every iteration. The per-row popcount of list_ep is not needed because
    // |size - hamming| is computed from the pair's hamming directly.
    std::vector<uint64_t> ep_bits(static_cast<size_t>(rows) * words, 0ULL);
    std::vector<int>      ep_pop(rows, 0); // unused; kept to share pack_bits()

    // Reusable ID buffers for the two halves of each bipartition.
    std::vector<int> ids_in, ids_out;
    ids_in.reserve(size);
    ids_out.reserve(size);

    #ifdef _OPENMP
    mt19937_64 mt(static_cast<uint64_t>(seed) + omp_get_thread_num() * 1000003ULL);
    #else
    mt19937_64 mt(static_cast<uint64_t>(seed));
    #endif
    uniform_real_distribution<double> urd(0, 1);

    // chunk=4: each iteration runs a full GS, so per-iteration work is very
    // large and roughly uniform — the default chunk=1 wastes scheduler
    // overhead on each handout.  4 keeps load-balance flexibility while
    // amortising the OpenMP runtime cost.
    #pragma omp for schedule(dynamic, 4)
    for(int n = 0; n < ep_num; n++){
      auto E_buf   = perturb(W, N, urd, mt);
      auto list_ep = run_gs_list(E_buf.get(), size);

      // Pack this iteration's list_ep into ep_bits (reused buffer).
      if(rows > 0){
        pack_bits(list_ep.get(), rows, size, words,
                  ep_bits.data(), ep_pop.data());
      }

      for(int x = 0; x < rows; x++){
        // Ori row's popcount is precomputed; minority count of branch x.
        int p_int = ori_pop[x];
        if(size - p_int < p_int) p_int = size - p_int;
        double p = static_cast<double>(p_int);

        // Find the minimum (size-corrected) Hamming distance to any ep row.
        // Early-exit at the theoretical minimum (delta == 0) since we cannot
        // improve below zero.
        int delta_int = size;
        const uint64_t* ori_row = ori_bits.data()
                                  + static_cast<size_t>(x) * words;
        for(int y = 0; y < rows; y++){
          const uint64_t* ep_row = ep_bits.data()
                                   + static_cast<size_t>(y) * words;
          int ham = 0;
          for(int w = 0; w < words; w++){
            ham += __builtin_popcountll(ori_row[w] ^ ep_row[w]);
          }
          if(size - ham < ham) ham = size - ham;
          if(ham < delta_int){
            delta_int = ham;
            if(delta_int == 0) break;
          }
        }
        double delta = static_cast<double>(delta_int);

        // Build the two clade keys via FNV-1a hash of sorted ID lists.
        if(delta < p - 1){
          ids_in.clear();
          ids_out.clear();
          int* ori_row_int = list_ori + static_cast<size_t>(x) * size;
          for(int z = 0; z < size; z++){
            if(ori_row_int[z] > 0) ids_in.push_back(z+1);
            else                   ids_out.push_back(z+1);
          }
          double const w = 1.0 - delta / (p - 1);
          local_ep[clade_key::hash_sorted_ids(ids_in)]  += w;
          local_ep[clade_key::hash_sorted_ids(ids_out)] += w;
        }
      }
    } // end omp for

    #pragma omp critical
    {
      for(auto const& kv : local_ep) ep[kv.first] += kv.second;
    }
  } // end parallel
}
