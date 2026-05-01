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
#include <memory>
#include <vector>
#include <cblas.h>
#include "transitivity.h"
#ifdef HAVE_CUDA
#  include "transitivity_gpu.h"
// Use GPU matrix multiply above this sequence count.
static constexpr int GPU_TRANS_THRESHOLD = 256;
#endif

// CPU sparse-path threshold.  When the off-diagonal density of the adjacency
// graph drops below this, the node-iterator triangle count is faster than the
// dense dsyrk path; above it, BLAS wins on raw throughput.  Empirically
// chosen — the crossover depends on machine BLAS performance, but 0.30 is a
// safe centre that avoids a regression on either side.
static constexpr double SPARSE_DENSITY_THRESHOLD = 0.30;

using namespace std;

// Sparse triangle-count formulation of the transitivity / global clustering
// coefficient.  Equivalent to the dense formula:
//     transitivity = Σ_{i<j} B[i,j]·A[i,j] / Σ_{i<j} B[i,j]
//                  = 3·#triangles / Σ_v C(deg(v), 2)
// (See e.g. Newman, "Networks", §7.9; the algebraic identities follow from
// B = A·A and A diagonal = 0.)
//
// Cost: O(Σ_v deg(v)²) for the triangle count + O(N²) to build the adjacency
// list (we have to scan W once anyway).  For a graph with average degree d,
// this is O(N·d²) — quadratic in d, not in N — so it dominates dsyrk by a
// factor of (N/d)·BLAS_efficiency on sparse inputs.
static double transitivity_sparse(double* const oW, int const N){
  // Build sorted adjacency lists in one pass.
  vector<vector<int>> adj(N);
  for(int i = 0; i < N; i++){
    auto& ai = adj[i];
    ai.reserve(16);
    double const* row = oW + static_cast<size_t>(i) * N;
    for(int j = 0; j < N; j++){
      if(j != i && row[j] > 0) ai.push_back(j); // ascending j ⇒ already sorted
    }
  }

  // Denominator: Σ_v C(deg(v), 2)
  long long denom = 0;
  for(int v = 0; v < N; v++){
    long long d = static_cast<long long>(adj[v].size());
    denom += d * (d - 1) / 2;
  }

  // Triangle count via node-iterator with strict ordering u < v < w.  Each
  // triangle {u,v,w} (u<v<w) is enumerated exactly once: when processing the
  // edge (u,v), we count common neighbours w > v.  Sorted-list intersection
  // restricted to the suffix > v.
  long long triangles = 0;
  for(int u = 0; u < N; u++){
    auto const& Nu = adj[u];
    for(int v : Nu){
      if(v <= u) continue;
      auto const& Nv = adj[v];
      auto it_u = std::lower_bound(Nu.begin(), Nu.end(), v + 1);
      auto it_v = std::lower_bound(Nv.begin(), Nv.end(), v + 1);
      while(it_u != Nu.end() && it_v != Nv.end()){
        if(*it_u == *it_v){
          triangles++;
          ++it_u; ++it_v;
        }
        else if(*it_u < *it_v) ++it_u;
        else                    ++it_v;
      }
    }
  }

  if(denom == 0) return 0.0;
  return 3.0 * static_cast<double>(triangles) / static_cast<double>(denom);
}

double transitivity(double* const (&oW), int const& size){

#ifdef HAVE_CUDA
  if (size > GPU_TRANS_THRESHOLD) {
    // Build adjacency matrix A (both triangles, row-major) then delegate to GPU
    auto A_buf = unique_ptr<double[]>(new double[size*size]());
    double* A  = A_buf.get();
    for(int i = 0; i < size-1; i++){
      for(int j = i+1; j < size; j++){
        double adj      = (oW[i*size+j] > 0) ? 1.0 : 0.0;
        A[i*size+j]     = adj;
        A[j*size+i]     = adj;
      }
    }
    return transitivity_gpu(A, size);
  }
#endif

  // Cheap density probe (one upper-triangle scan, O(N²/2)).  Cost is dwarfed
  // by either downstream path, so it's worth paying to dispatch correctly.
  long long edges = 0;
  for(int i = 0; i < size-1; i++){
    double const* row = oW + static_cast<size_t>(i) * size;
    for(int j = i+1; j < size; j++){
      if(row[j] > 0) edges++;
    }
  }
  long long const max_pairs = static_cast<long long>(size) * (size - 1) / 2;
  double density = (max_pairs > 0)
                 ? static_cast<double>(edges) / static_cast<double>(max_pairs)
                 : 0.0;

  if(density < SPARSE_DENSITY_THRESHOLD){
    return transitivity_sparse(oW, size);
  }

  // CPU dense path (cblas_dsyrk) ────────────────────────────────────────
  auto A_buf = unique_ptr<double[]>(new double[size*size]());
  auto B_buf = unique_ptr<double[]>(new double[size*size]());
  double* A = A_buf.get();
  double* B = B_buf.get();

  // A: adjacency matrix of oW (both triangles)
  for(int i=0; i<size-1; i++){
    for(int j=i+1; j<size; j++){
      double adj = (oW[i*size+j]>0)? 1: 0;
      A[i*size+j] = adj;
      A[i+size*j] = adj;
    }
  }

  // B = A * A^T  (upper triangle via cblas_dsyrk)
  const double ALPHA = 1;
  const double BETA  = 0;
  cblas_dsyrk(
      CblasRowMajor,
      CblasUpper,
      CblasNoTrans,
      size, size,
      ALPHA, A, size,
      BETA,  B, size
  );

  // transitivity = Σ B[i][j]*A[i][j] / Σ B[i][j]  (i < j)
  double numerator   = 0;
  double denominator = 0;
  for(int i=0; i<size-1; i++){
    for(int j=i+1; j<size; j++){
      numerator   += B[i*size+j]*A[i*size+j];
      denominator += B[i*size+j];
    }
  }

  if(denominator > 0){
    return numerator/denominator;
  }
  return 0;
}
