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

#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <iostream>
#include <memory>
#include "transitivity_gpu.h"

using namespace std;

static cublasHandle_t tl_cublas = nullptr;

static cublasHandle_t get_cublas() {
    if (!tl_cublas) {
        if (cublasCreate(&tl_cublas) != CUBLAS_STATUS_SUCCESS) {
            tl_cublas = nullptr;
        }
    }
    return tl_cublas;
}

// -----------------------------------------------------------------------
// Layout note (adjacency matrix A, row-major, both triangles filled)
// -----------------------------------------------------------------------
// cuBLAS uses column-major storage.  We want B = A * A^T.
// Since A is symmetric (A = A^T), B = A * A = A^2.
//
// cublasDsyrk(FILL_MODE_LOWER, OP_N) computes the column-major lower
// triangle of  C = alpha * A * A^T + beta * C.
// Column-major lower triangle ↔ row-major upper triangle.
// ⇒ After copying d_B to host, B_host[i*N + j] for j > i holds B[i][j]. ✓
// -----------------------------------------------------------------------

double transitivity_gpu(double* const A, int const size) {
    cublasHandle_t handle = get_cublas();
    if (!handle) {
        cerr << "[transitivity_gpu] cuBLAS handle unavailable" << endl;
        return 0.0;
    }

    double* d_A = nullptr;
    double* d_B = nullptr;
    cudaMalloc(&d_A, sizeof(double) * size * size);
    cudaMalloc(&d_B, sizeof(double) * size * size);
    cudaMemset(d_B, 0, sizeof(double) * size * size);

    cudaMemcpy(d_A, A, sizeof(double) * size * size, cudaMemcpyHostToDevice);

    // B = 1 * A * A^T + 0 * B  (fills column-major lower = row-major upper)
    const double alpha = 1.0, beta = 0.0;
    cublasDsyrk(handle,
                CUBLAS_FILL_MODE_LOWER, // col-major lower = row-major upper
                CUBLAS_OP_N,
                size, size,
                &alpha, d_A, size,
                &beta,  d_B, size);

    auto B_host = unique_ptr<double[]>(new double[size * size]());
    cudaMemcpy(B_host.get(), d_B, sizeof(double) * size * size,
               cudaMemcpyDeviceToHost);

    cudaFree(d_A);
    cudaFree(d_B);

    // transitivity = Σ B[i][j]*A[i][j]  /  Σ B[i][j]   (i < j only)
    double numerator   = 0.0;
    double denominator = 0.0;
    for (int i = 0; i < size - 1; i++) {
        for (int j = i + 1; j < size; j++) {
            double b = B_host[i * size + j]; // row-major upper triangle ✓
            numerator   += b * A[i * size + j];
            denominator += b;
        }
    }

    return (denominator > 0.0) ? numerator / denominator : 0.0;
}
