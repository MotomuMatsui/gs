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

tuple<vector<int>,vector<int>> spectral_clustering(double* const& oW, vector<int> const& res, int const& num){

  // Cluster size
  int N = (int)count(res.begin(), res.end(), num);

  // Owned buffers; raw aliases allow passing to legacy double*(&) interfaces.
  // W and A are fully written below before being read, so skip zero-init.
  // D is an accumulator (D[k] += w) so must start at zero.
  auto W_buf = unique_ptr<double[]>(new double[N*N]);
  auto D_buf = unique_ptr<double[]>(new double[N]());
  auto A_buf = unique_ptr<double[]>(new double[N*N]);
  double* W = W_buf.get();
  double* D = D_buf.get();
  double* A = A_buf.get();

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
  // write exactly one length-N column.  We allocate just that column — N
  // doubles instead of N*N — saving size² alloc churn per cluster split.
  auto z_buf = unique_ptr<double[]>(new double[N]); // Fiedler vector buffer
  double* z  = z_buf.get();
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

  //Split submatrix into two parts; whichCUT writes all N entries of qi.
  auto qi_buf = unique_ptr<int[]>(new int[N]);
  int* qi = qi_buf.get();
  int  cut = 0;
  vector<int> a, b;
  whichCUT(z, D, qi, cut, N);
  splitVECTOR(res, num, qi, cut, a, b, N);

  return make_tuple(std::move(a), std::move(b));
}
