/******************************************\
| Graph Splitting Method v2.2 (2018/11/07) |
|                                          |
| Copyright (c) 2015-2018 Motomu Matsui    |
|     Distributed under the GNU GPL        |
|                                          |
|     Matsui M and Iwasaki W (2018)        |
|     Systematic Biology, xx:xx-xx.        |
|                                          |
|     http://gs.bs.s.u-tokyo.ac.jp/        |
\******************************************/

#include <lapacke.h>

using namespace std;

int eigen_lapack(double* (&A), double* (&z), int N){

  //Variables for LAPACKE_dsyevr function
  double  vl     = 1;
  double  vu     = N;
  int     il     = 2; // Graph Splitting method (Spectral clustering)
  int     iu     = 2; // requires only the 2nd smallest eigenvalue
  double  abstol = 0;

  int* m      = new int[iu-il+1];
  double* w   = new double[N];
  int* isuppz = new int[(iu-il+1)*N];
  
  //Eigenvalue analysis using LAPACK/BLAS (real symmetric matrix)
  auto info = LAPACKE_dsyevr(
			     LAPACK_ROW_MAJOR, // Matrix layout (Row major layout)
			     'V',              // Jobs (V: eigenvalues and eigenvectors)
			     'I',              // Range (I: the IL-th through IU-th eigenvalues will be found)
			     'L',              // Stored triangle matrix (L: Lower triangle of A is stored)
			     N,                // Order of A
			     A,                // Real symmetric matrix
			     N,                // Leading dimension of A
			     vl,vu,            // (These arguments are NOT referenced, but required as dummies)
			     il,iu,            // Index of the smallest (il) and largest (iu) eigenvalues to be returned
			     abstol,           // Absolute error tolerance
			     m,                // Total number of eigenvalues
			     w,                // Eivenvalues
			     z,                // Eigenvectors
			     N,                // Leading dimension of z
			     isuppz            // Support value of z
			     );
  
  if(info != 0){
    //Error occurs if info>0
  }

  //Free
  delete[] m;
  delete[] w;
  delete[] isuppz;

  return info;
}
