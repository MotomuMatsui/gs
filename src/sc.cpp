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
#include <cstring>
#include <iostream>
#include <memory>
#include <tuple>
#include <vector>

#include "eigen.h"
#ifdef HAVE_CUDA
#  include "eigen_gpu.h"
// Use GPU eigenvalue solver when the submatrix exceeds this size.
// Below the threshold, PCIe transfer overhead dominates.
static constexpr int GPU_EIGEN_THRESHOLD = 64;
#endif

#include "sc.h"
#include "sc_functions.h"

using namespace std;

namespace {
// Per-thread scratch reused across cluster splits.  Each thread (the main
// thread for plain GS, every OMP worker for EP) grows its buffers up to the
// largest cluster size it has ever seen, then reuses them — eliminating the
// `new double[N*N]` × 5 churn that bottlenecked EP at high OMP thread counts
// (each split paying glibc-arena and page-fault costs).
struct SCScratch {
  vector<double> W;   // N*N
  vector<double> A;   // N*N
  vector<double> D;   // N (accumulator — must be zeroed before subMATRIX)
  vector<double> z;   // N (Fiedler vector)
  vector<int>    qi;  // N (permutation written by whichCUT)
  int capacity = 0;

  void reserve(int N) {
    if (N <= capacity) return;
    // Grow exponentially-ish to amortise the resize cost across splits within
    // one EP iter (cluster sizes shrink monotonically anyway, so we usually
    // grow only on the very first call).
    size_t sq = static_cast<size_t>(N) * N;
    W.resize(sq);
    A.resize(sq);
    D.resize(N);
    z.resize(N);
    qi.resize(N);
    capacity = N;
  }
};
thread_local SCScratch sc_scratch;
} // namespace

tuple<vector<int>,vector<int>> spectral_clustering(double* const& oW, vector<int> const& res, int const& num){

  // Cluster size
  int N = (int)count(res.begin(), res.end(), num);

  sc_scratch.reserve(N);
  double* W = sc_scratch.W.data();
  double* A = sc_scratch.A.data();
  double* D = sc_scratch.D.data();
  double* z = sc_scratch.z.data();
  int*    qi = sc_scratch.qi.data();
  // D is an accumulator (subMATRIX does D[k] += w); zero the first N entries
  // each call.  W and A are fully overwritten below before being read.
  std::memset(D, 0, sizeof(double) * static_cast<size_t>(N));

  //Rayleigh quotient of the submatrix
  subMATRIX(oW, W, D, res, num);

  //A = E - D*W*D  (normalised Laplacian, lower triangle)
  for(int x = 0; x < N; x++){
    for(int y = 0; y < x; y++){
      auto p  = x*N+y;
      A[p]    = -D[x]*D[y]*W[p];
    }
    auto p  = x*N+x;
    A[p]    = 1-D[x]*D[x]*W[p];
  }

  //Eigenvalue analysis: choose CPU or GPU path.
  // Only the Fiedler vector is consumed by whichCUT, and both backends
  // (eigen_lapack with il=iu=2, eigen_gpu reading column 1 of cuSOLVER output)
  // write exactly one length-N column.  `z` is the thread-local Fiedler
  // buffer (N doubles) drawn from sc_scratch above.
  int info   = 0;

#ifdef HAVE_CUDA
  if (N > GPU_EIGEN_THRESHOLD) {
    info = eigen_gpu(A, z, N);
    if (info != 0) {
      // GPU failed — fall back to CPU LAPACK.  Zero-fill is defensive; both
      // backends overwrite the entire buffer on success.
      fill(z, z + N, 0.0);
      info = eigen_lapack(A, z, N);
    }
  } else {
    info = eigen_lapack(A, z, N);
  }
#else
  info = eigen_lapack(A, z, N);
#endif

  if(info != 0){
    // On unrecoverable failure: return entire cluster as 'a', empty 'b'
    vector<int> a, b;
    for(int i = 0; i < (int)res.size(); i++){
      if(res[i] == num) a.push_back(i);
    }
    return make_tuple(std::move(a), std::move(b));
  }

  //Split submatrix into two parts; whichCUT writes all N entries of qi
  //into the thread-local scratch buffer above.
  int  cut = 0;
  vector<int> a, b;
  whichCUT(z, D, qi, cut, N);
  splitVECTOR(res, num, qi, cut, a, b, N);

  return make_tuple(std::move(a), std::move(b));
}
