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

#ifndef TRANSITIVITY_GPU_H
#define TRANSITIVITY_GPU_H

// GPU-accelerated transitivity using cuBLAS dsyrk.
// A : size×size adjacency matrix (row-major, both triangles filled)
double transitivity_gpu(double* const A, int const size);

#endif /*TRANSITIVITY_GPU_H*/
