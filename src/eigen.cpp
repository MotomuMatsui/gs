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

#include <iostream>
#include <memory>
#include <lapacke.h>
#include "eigen.h"

using namespace std;

int eigen_lapack(double* (&A), double* (&z), int N){

  //Variables for LAPACKE_dsyevr function
  double vl     = 1;
  double vu     = N;
  int    il     = 2; // Graph Splitting method (Spectral clustering)
  int    iu     = 2; // requires only the 2nd smallest eigenvalue (Fiedler vector)
  double abstol = 0; // absolute error tolerance (0 = use default machine precision)

  auto m      = unique_ptr<int[]>   (new int   [iu-il+1]());
  auto w      = unique_ptr<double[]>(new double [N]());
  auto isuppz = unique_ptr<int[]>   (new int   [(iu-il+1)*N]());

  //Eigenvalue analysis using LAPACK/BLAS (real symmetric matrix)
  // ldz=1: only one eigenvector (M=iu-il+1=1) is requested, and we store it
  // contiguously in row-major as z[0..N-1].  Passing ldz=N would force LAPACKE
  // to expect an N×N output buffer; we only have N doubles allocated.
  auto info = LAPACKE_dsyevr(
      LAPACK_ROW_MAJOR, // Matrix layout (Row major layout)
      'V',              // Jobs (V: eigenvalues and eigenvectors)
      'I',              // Range (I: the IL-th through IU-th eigenvalues will be found)
      'L',              // Stored triangle matrix (L: Lower triangle of A is stored)
      N,                // Order of A
      A,                // Real symmetric matrix
      N,                // Leading dimension of A
      vl, vu,           // (These arguments are NOT referenced when range='I', but required as dummies)
      il, iu,           // Index of the smallest (il) and largest (iu) eigenvalues to be returned
      abstol,           // Absolute error tolerance
      m.get(),          // Total number of eigenvalues found
      w.get(),          // Eigenvalues
      z,                // Eigenvectors (length N, one column)
      1,                // Leading dimension of z (row stride in row-major; M=1)
      isuppz.get()      // Support value of z
  );

  if(info != 0){
    cerr << "[eigen_lapack] LAPACK dsyevr failed: info=" << info
         << " (N=" << N << ")" << endl;
  }

  return info;
}
