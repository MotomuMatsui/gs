/******************************************\
| Graph Splitting Method v1.0 (2018/06/11) |
|                                          |
| Copyright (c) 2015-2018 Motomu Matsui    |
|     Distributed under the GNU GPL        |
|                                          |
|     Matsui M and Iwasaki W (2018)        |
|     Systematic Biology, xx:xx-xx.        |
|                                          |
|     http://gs.bs.s.u-tokyo.ac.jp/        |
\******************************************/

#include <vector>
#include <tuple>

using namespace std;

//sc.cpp (spectral clustering)
extern tuple<vector<int>,vector<int>> spectral_clustering(double* const&, vector<int> const&, int const&);

//gs_functions.cpp
extern void sedMATRIX(int*&, vector<int>&, vector<int> const&, int const, int const);
extern double simI(double* const&, vector<int> const&, int const);
extern int whichMIN(vector<double> const&);

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
