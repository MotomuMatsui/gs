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

#ifndef EIGEN_GPU_H
#define EIGEN_GPU_H

// True if at least one CUDA device is present
bool cuda_available();

// GPU Fiedler-vector solver (cuSOLVER dsyevd).
// A : N×N, row-major, lower triangle filled (symmetric Laplacian)
// z : N×N output — Fiedler vector written into column 0 (stride N),
//     matching the layout expected by whichCUT(z, col=0, ...).
// Returns the cuSOLVER info code (0 = success).
int eigen_gpu(double* A, double* z, int N);

#endif /*EIGEN_GPU_H*/
