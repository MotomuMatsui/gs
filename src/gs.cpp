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

#include <memory>
#include <tuple>
#include <vector>

#include "gs.h"
#include "gs_functions.h"
#include "sc.h"

using namespace std;

void GS(double* const (&W), int* (&step), int const& size){

  vector<double> simL(size, 3); // 3 > 1 (max similarity), acts as sentinel for unvisited clusters
  vector<int> res(size, 1);

  // Manage step internally; ownership transfers to caller via release()
  auto step_buf = unique_ptr<int[]>(new int[size*size]());
  int* step_raw = step_buf.get();
  sedMATRIX(step_raw, res, res, 1, size);

  int gK   = 1;
  int gMin = 1;
  vector<int> a;
  vector<int> b;

  while(gK < size){
    //Spectral clustering
    tie(a,b) = spectral_clustering(W, res, gMin);

    //Record the result of spectral clustering
    gK++;
    sedMATRIX(step_raw, res, b, gK, size);

    //Decision of the cluster which will be analyzed at the next step
    simL[gMin-1] = simI(W, res, gMin); // cluster a
    simL[gK-1]   = simI(W, res, gK);   // cluster b

    gMin = whichMIN(simL);
  }

  // Transfer ownership to caller
  step = step_buf.release();
}
