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

#include <cuda_runtime.h>
#include <cusolverDn.h>
#include <iostream>
#include <vector>
#include "eigen_gpu.h"

using namespace std;

// One handle per thread so OpenMP threads don't share state.
// Handles are lazily created and never explicitly destroyed
// (CUDA runtime cleans up on process exit).
static thread_local cusolverDnHandle_t tl_solver = nullptr;

static cusolverDnHandle_t get_solver() {
    if (!tl_solver) {
        if (cusolverDnCreate(&tl_solver) != CUSOLVER_STATUS_SUCCESS) {
            tl_solver = nullptr;
        }
    }
    return tl_solver;
}

bool cuda_available() {
    int n = 0;
    return (cudaGetDeviceCount(&n) == cudaSuccess) && (n > 0);
}

// -----------------------------------------------------------------------
// Layout note (symmetric matrix A, row-major lower triangle filled)
// -----------------------------------------------------------------------
// cuSOLVER uses column-major storage.  For a symmetric A stored row-major
// the flat-array contents are identical to the transposed matrix in column-
// major, i.e. cuSOLVER sees A^T.  Since A is symmetric, A^T = A, so the
// eigenvalues and eigenvectors are unchanged.
//
// Which fill mode to use?
//   Row-major lower triangle  ↔  column-major upper triangle.
//   (flat index r*N+c with r>=c equals column-major index r+c*N only when
//    we map r→col, c→row, giving col>=row, i.e. upper in col-major.)
// ⇒ CUBLAS_FILL_MODE_UPPER reads exactly the filled entries of A.
// -----------------------------------------------------------------------

// Helper: free any non-null device pointer in argument list.
static void free_all(double* a, double* b, double* c, int* d) {
    if (a) cudaFree(a);
    if (b) cudaFree(b);
    if (c) cudaFree(c);
    if (d) cudaFree(d);
}

int eigen_gpu(double* A, double* z, int N) {
    cusolverDnHandle_t handle = get_solver();
    if (!handle) {
        cerr << "[eigen_gpu] cuSOLVER handle unavailable" << endl;
        return -1;
    }

    double* d_A    = nullptr;
    double* d_W    = nullptr;
    double* d_work = nullptr;
    int*    d_info = nullptr;
    int     lwork  = 0;

    cudaError_t cerr_code;

    cerr_code = cudaMalloc(&d_A, sizeof(double) * N * N);
    if (cerr_code != cudaSuccess) {
        cerr << "[eigen_gpu] cudaMalloc(d_A) failed: "
             << cudaGetErrorString(cerr_code) << " (N=" << N << ")" << endl;
        free_all(d_A, d_W, d_work, d_info);
        return -1;
    }
    cerr_code = cudaMalloc(&d_W, sizeof(double) * N);
    if (cerr_code != cudaSuccess) {
        cerr << "[eigen_gpu] cudaMalloc(d_W) failed: "
             << cudaGetErrorString(cerr_code) << endl;
        free_all(d_A, d_W, d_work, d_info);
        return -1;
    }
    cerr_code = cudaMalloc(&d_info, sizeof(int));
    if (cerr_code != cudaSuccess) {
        cerr << "[eigen_gpu] cudaMalloc(d_info) failed: "
             << cudaGetErrorString(cerr_code) << endl;
        free_all(d_A, d_W, d_work, d_info);
        return -1;
    }

    // Upload row-major A; cuSOLVER will interpret it as column-major A^T = A
    cerr_code = cudaMemcpy(d_A, A, sizeof(double) * N * N, cudaMemcpyHostToDevice);
    if (cerr_code != cudaSuccess) {
        cerr << "[eigen_gpu] cudaMemcpy(H2D) failed: "
             << cudaGetErrorString(cerr_code) << endl;
        free_all(d_A, d_W, d_work, d_info);
        return -1;
    }

    // Query workspace
    cusolverStatus_t st = cusolverDnDsyevd_bufferSize(
        handle,
        CUSOLVER_EIG_MODE_VECTOR,
        CUBLAS_FILL_MODE_UPPER, // reads row-major lower tri
        N, d_A, N, d_W, &lwork);
    if (st != CUSOLVER_STATUS_SUCCESS) {
        cerr << "[eigen_gpu] cusolverDnDsyevd_bufferSize failed: status=" << st << endl;
        free_all(d_A, d_W, d_work, d_info);
        return -1;
    }

    cerr_code = cudaMalloc(&d_work, sizeof(double) * lwork);
    if (cerr_code != cudaSuccess) {
        cerr << "[eigen_gpu] cudaMalloc(d_work) failed: "
             << cudaGetErrorString(cerr_code) << " (lwork=" << lwork << ")" << endl;
        free_all(d_A, d_W, d_work, d_info);
        return -1;
    }

    // Eigendecomposition — eigenvectors overwrite d_A (column-major, sorted ASC)
    st = cusolverDnDsyevd(handle,
                          CUSOLVER_EIG_MODE_VECTOR,
                          CUBLAS_FILL_MODE_UPPER,
                          N, d_A, N, d_W, d_work, lwork, d_info);
    if (st != CUSOLVER_STATUS_SUCCESS) {
        cerr << "[eigen_gpu] cusolverDnDsyevd failed: status=" << st << endl;
        free_all(d_A, d_W, d_work, d_info);
        return -1;
    }

    int h_info = 0;
    cerr_code = cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost);
    if (cerr_code != cudaSuccess) {
        cerr << "[eigen_gpu] cudaMemcpy(info D2H) failed: "
             << cudaGetErrorString(cerr_code) << endl;
        free_all(d_A, d_W, d_work, d_info);
        return -1;
    }

    if (h_info == 0) {
        // Fiedler vector = 2nd smallest eigenvector = column 1 (0-indexed)
        // in column-major d_A → contiguous at d_A + N, length N.
        // whichCUT expects z[p*N] for p=0..N-1 (col=0, stride=N).
        vector<double> tmp(N);
        cerr_code = cudaMemcpy(tmp.data(), d_A + N, sizeof(double) * N,
                               cudaMemcpyDeviceToHost);
        if (cerr_code != cudaSuccess) {
            cerr << "[eigen_gpu] cudaMemcpy(fiedler D2H) failed: "
                 << cudaGetErrorString(cerr_code) << endl;
            free_all(d_A, d_W, d_work, d_info);
            return -1;
        }
        for (int p = 0; p < N; p++) z[p * N] = tmp[p];
    } else {
        cerr << "[eigen_gpu] cusolverDnDsyevd failed: info=" << h_info
             << " (N=" << N << ")" << endl;
    }

    free_all(d_A, d_W, d_work, d_info);
    return h_info;
}
