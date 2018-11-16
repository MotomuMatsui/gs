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
extern void sedMATRIX(int*&, vector<int>&, vector<int> const&, int const, int const);
extern double simI(double* const&, vector<int> const&, int const);
extern int whichMIN(vector<double> const&);
extern double gev(double const&, double const&);

//format.cpp
extern void sc2list(int* const&, int*&, int const&);

void EP(double* const (&W), unordered_map<string, double>& ep, function<double()>& R, int const& size){

  // Edge perturbation
  auto N = size*size;
  auto E = new double[N]();
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

  auto step = new int[N]();
  sedMATRIX(step, res, res, 1, size);

  int gK   = 1;
  int gMin = 1;
  vector<int> a;
  vector<int> b;

  while(gK < size){
    //Spectral clustering
    tie(a,b) = spectral_clustering(E, res, gMin);

    //Record the result of spectral clustering
    gK++;
    sedMATRIX(step, res, b, gK, size);
    
    //Decision of the cluster which will be analyzed at the next step
    simL[gMin-1] = simI(E, res, gMin); //cluster a
    simL[gK-1]   = simI(E, res, gK);   //cluster b

    gMin = whichMIN(simL);
  }

  int* list_ep;
  sc2list(step, list_ep, size);

  stringstream ss;
  for(int x=0; x<2*(size-3); x++){
    ss.str("");
    for(int z=0; z<size; z++){
      if(list_ep[x*size+z]>0){
        ss << z+1 << "|";
      }
    }
    ep[ss.str()] ++;
  }

  //Free
  delete[] step;
  delete[] list_ep;
  delete[] E;
}

void EP2(double* const (&W), int* const (&list), unordered_map<string, double>& ep, function<double()>& R, int const& size){

  // Edge perturbation
  auto N = size*size;
  auto E = new double[N]();
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

  auto step = new int[N]();
  sedMATRIX(step, res, res, 1, size);

  int gK   = 1;
  int gMin = 1;
  vector<int> a;
  vector<int> b;

  while(gK < size){
    //Spectral clustering
    tie(a,b) = spectral_clustering(E, res, gMin);

    //Record the result of spectral clustering
    gK++;
    sedMATRIX(step, res, b, gK, size);

    //Decision of the cluster which will be analyzed at the next step
    simL[gMin-1] = simI(E, res, gMin); //cluster a
    simL[gK-1]   = simI(E, res, gK);   //cluster b

    gMin = whichMIN(simL);
  }

  int* list_ep;
  sc2list(step, list_ep, size);

  for(int x=0; x<2*(size-3); x++){
    double min_share = size;
    for(int y=0; y<2*(size-3); y++){
      double share = 0;
      int num_a = 0;
      int num_b = 0;
      for(int z=0; z<size; z++){
        share += list[x*size+z]^list_ep[y*size+z];
	num_a += list[x*size+z];
	num_b += list_ep[y*size+z];
      }
      if(num_a == num_b){
	min_share = (share < min_share)? share: min_share;
      }
    }
    
    stringstream ss;
    int num = 0;
    for(int z=0; z<size; z++){
      if(list[x*size+z]>0){
        ss << z+1 << "|";       
        num ++;
      }
    }
    if(num > 1 && min_share<num-1){
      ep[ss.str()] += 1 - min_share / (num - 1);
    }
  }

  //Free
  delete[] step;
  delete[] list_ep;
  delete[] E;
}
