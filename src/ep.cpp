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
#include <cmath> 
#include <functional>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#include <unordered_map>

#include "sc.h"
#include "gs_functions.h"
#include "format.h"

using namespace std;

void EP_fbs(double* const (&W), unordered_map<string, double>& ep, function<double()>& R, int const& size){

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

  stringstream ss1;
  stringstream ss2;
  for(int x=0; x<size-3; x++){
    ss1.str("");
    ss2.str("");
    for(int z=0; z<size; z++){
      if(list_ep[x*size+z]>0){
        ss1 << z+1 << "|";
      }
      else{
	ss2 << z+1 << "|";
      }
    }
    ep[ss1.str()] ++;
    ep[ss2.str()] ++;
  }

  //Free
  delete[] step;
  delete[] list_ep;
  delete[] E;
}

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

void EP_tbe(double* const (&W), int* const (&list_ori), unordered_map<string, double>& ep, function<double()>& R, int const& size){

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

  for(int x=0; x<size-3; x++){
    // calculate subset sizes
    double p = 0;
    for(int z=0; z<size; z++){
      p += list_ori[x*size+z];
    }
    p = (p <= size - p)? p: size - p;

    double delta = size;
    for(int y=0; y<size-3; y++){
      double hamming = 0;
      for(int z=0; z<size; z++){
	hamming += list_ori[x*size+z]^list_ep[y*size+z];
      }
      hamming = (hamming < size - hamming)? hamming: size - hamming;

      delta = (hamming < delta)? hamming: delta;
    }

    stringstream ss1;
    stringstream ss2;
    for(int z=0; z<size; z++){
      if(list_ori[x*size+z]>0){
        ss1 << z+1 << "|";      
      }
      else{
        ss2 << z+1 << "|";      
      }
    }
    if(delta < p - 1){
      ep[ss1.str()] += 1 - delta / (p - 1);
      ep[ss2.str()] += 1 - delta / (p - 1);
    }
  }

  //Free
  delete[] step;
  delete[] list_ep;
  delete[] E;
}

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
