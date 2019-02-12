/******************************************\
| Graph Splitting Method v2.4 (2019/02/12) |
|                                          |
| Copyright (c) 2015-2019 Motomu Matsui    |
|     Distributed under the GNU GPL        |
|                                          |
|     Matsui M and Iwasaki W (2019)        |
|     Systematic Biology, xx:xx-xx.        |
|                                          |
|     http://gs.bs.s.u-tokyo.ac.jp/        |
\******************************************/

#include <vector>
#include <tuple>

#include "gs.h"
#include "gs_functions.h"
#include "sc.h"

using namespace std;

void GS(double* const (&W), int* (&step), int const& size){

  vector<double> simL(size,3);
  vector<int> res(size,1);

  step = new int[size*size]();
  sedMATRIX(step, res, res, 1, size);

  int gK   = 1;
  int gMin = 1;
  vector<int> a;
  vector<int> b;

  while(gK < size){
    //Spectral clustering
    tie(a,b) = spectral_clustering(W, res, gMin);

    //Record the result of spectral clustering
    gK++;
    sedMATRIX(step, res, b, gK, size);
    
    //Decision of the cluster which will be analyzed at the next step
    simL[gMin-1] = simI(W, res, gMin); //cluster a
    simL[gK-1]   = simI(W, res, gK);   //cluster b

    gMin = whichMIN(simL);
  }
}
