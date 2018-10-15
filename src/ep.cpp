/******************************************\
| Graph Splitting Method v2.1 (2018/10/15) |
|                                          |
| Copyright (c) 2015-2018 Motomu Matsui    |
|     Distributed under the GNU GPL        |
|                                          |
|     Matsui M and Iwasaki W (2018)        |
|     Systematic Biology, xx:xx-xx.        |
|                                          |
|     http://gs.bs.s.u-tokyo.ac.jp/        |
\******************************************/

#include <random>   
#include <cmath> 
#include <functional>
#include <algorithm>
#include <vector>
#include <tuple>
#include <sstream>
#include <string>
#include <unordered_map>
#include <iostream>

using namespace std;

//sc.cpp (spectral clustering)
extern tuple<vector<int>,vector<int>> spectral_clustering(double* const&, vector<int> const&, int const&);

//gs_functions.cpp
extern void sedVECTOR(vector<int>&, vector<int> const&, int const);
extern double simI(double* const&, vector<int> const&, int const);
extern int whichMIN(vector<double> const&);
extern double gev(double const&, double const&);

void EP(double* const (&W), unordered_map<string, double>& ep, function<double()>& R, int const& size){

  // Edge perturbation
  int N = size*size;
  double* E = new double[N]();
  for(int n = 0; n < N; n++){
    double a = W[n];
    if(a==0){
      E[n] = 0; // Graph topology is not changed
    }
    else{
      auto b = gev(R(), a);
      E[n] = 
	(b>1)? 1: // b: Similarity scores perturved according to Generalized Extreme Value distribution
	(b<0)? 0: // E[n]: Perturved sequence similarity score
	       b; // 0 <= E[n] <= 1
    }
  }

  //GS method + EP method
  vector<double> simL(size,3);
  vector<int> res(size,1);

  int gK   = 1;
  int gMin = 1;
  vector<int> a;
  vector<int> b;

  stringstream ss;

  while(gK < size){
    //Spectral clustering
    tie(a,b) = spectral_clustering(E, res, gMin);
    
    ss.str("");
    sort(a.begin(), a.end());
    for(int n : a){
      ss << n+1 << "|";
    }
    ep[ss.str()]++;

    ss.str("");
    sort(b.begin(), b.end());
    for(int n : b){
      ss << n+1 << "|";
    }
    ep[ss.str()]++;

    //Record the result of spectral clustering
    gK++;
    sedVECTOR(res, b, gK);
    
    //Decision of the cluster which will be analyzed at the next step
    simL[gMin-1] = simI(E, res, gMin); //cluster a
    simL[gK-1]   = simI(E, res, gK);   //cluster b

    gMin = whichMIN(simL);
  }

  //Free
  delete[] E;
}
