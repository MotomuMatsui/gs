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

#include <algorithm>
#include <tuple>
#include <vector>

#include "eigen.h"
#include "sc.h"
#include "sc_functions.h"

using namespace std;

tuple<vector<int>,vector<int>> spectral_clustering(double* const& oW, vector<int> const& res, int const& num){

  // Cluster size
  int N = (int)count(res.begin(), res.end(), num);
  
  double* W = new double[N*N]();
  double* D = new double[N]();
  double* A = new double[N*N]();  

  //Rayleigh quotient of the submatrix
  subMATRIX(oW, W, D, res, num);

  //A = E-D*W*D
  for(int x = 0; x < N; x++){
    for(int y = 0; y < x; y++){
      auto p = x*N+y;
      A[p] = -D[x]*D[y]*W[p];
    }
    auto p = x*N+x;
    A[p] = 1-D[x]*D[x]*W[p];
  }

  //Eigenvalue analysis
  double* z = new double[N*N]();     //z: eigenvector
  auto info = eigen_lapack(A, z, N); //LAPACK!
  if(info != 0){}                    //Error occurs if info>0

  //Split submatrix into two parts according to the result of spectral clustering
  int* qi = new int[N]();
  int cut = 0;
  vector<int> a, b;
  whichCUT(z, 0, D, qi, cut, N);
  splitVECTOR(res, num, qi, cut, a, b, N);
  
  //free
  delete[] A;
  delete[] W;
  delete[] D;
  delete[] z;
  delete[] qi;

  return forward_as_tuple(a,b);
}
