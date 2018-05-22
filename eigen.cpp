/*=========================================*\
| Graph Splitting Method v1.0 (2018/06/11) |
|                                          |
| Copyright (c) 2018 Motomu Matsui         |
| Systematic Biology, xx:xx-xx, 2018       |
|                                          |
| Web:  http://gs.bs.s.u-tokyo.ac.jp/      |
| Mail: qm at bs.s.u-tokyo.ac.jp           |
\*=========================================*/

#include <lapacke.h>

using namespace std;

int eigen_lapack(double* (&A), double* (&z), int N){

  //Variables for LAPACKE_dsyevr function
  int     info   = 0;
  double  vl     = 1;
  double  vu     = N;
  int     il     = 2; // Graph Splitting method (Spectral clustering)
  int     iu     = 2; // requires only the 2nd smallest eigenvalue
  double  abstol = 0;

  int* m      = new int[iu-il+1];
  double* w   = new double[N];
  int* isuppz = new int[(iu-il+1)*N];
  
  //Eigenvalue analysis using LAPACK/BLAS (real symmetric matrix)
  info = LAPACKE_dsyevr(
			LAPACK_ROW_MAJOR, // Matrix layout (Row major layout)
			'V',              // Jobs (V: eigenvalues and eigenvectors)
			'I',              // Range (I: the IL-th through IU-th eigenvalues will be found)
			'L',              // Stored triangle matrix (L: Lower triangle of A is stored)
			N,                // Order of A
			A,                // Real symmetric matrix
			N,                // Leading dimension of A
			vl,vu,            // (Not referenced, but this function requires dummy numbers)
			il,iu,            // Index of the smallest (il) and largest (iu) eigenvalue to be returned
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
