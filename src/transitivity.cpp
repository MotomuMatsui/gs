/******************************************\
| Graph Splitting Method v2.3 (2018/11/16) |
|                                          |
| Copyright (c) 2015-2018 Motomu Matsui    |
|     Distributed under the GNU GPL        |
|                                          |
|     Matsui M and Iwasaki W (2018)        |
|     Systematic Biology, xx:xx-xx.        |
|                                          |
|     http://gs.bs.s.u-tokyo.ac.jp/        |
\******************************************/

#include <cblas.h>
#include "transitivity.h"

using namespace std;

double transitivity(double* const (&oW), int const& size){

  double* A = new double[size*size]();
  double* B = new double[size*size]();

  // A: adjacency matrix of oW
  for(int i=0; i<size-1; i++){
    for(int j=i+1; j<size; j++){
      double adj = (oW[i*size+j]>0)? 1: 0;
      A[i*size+j] = adj;
      A[i+size*j] = adj;
    }
  }
  
  // B = 1*AA + 0*B
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

  // transitivity = sum(A^2*A)/sum(A^2)
  double numerator    = 0; //i-k&k-j&j-i (# of complete triangles)
  double denominator  = 0; //i-k&k-j&j?i (# of pathways whose lengths are 2)
  double transitivity = 0;
  
  for(int i=0; i<size-1; i++){
    for(int j=i+1; j<size; j++){
      numerator   += B[i*size+j]*A[i*size+j];
      denominator += B[i*size+j];
    }
  }
  
  if(denominator > 0){
    transitivity = numerator/denominator;
  }

  //Free
  delete[] A;
  delete[] B;
  
  return transitivity;
}
