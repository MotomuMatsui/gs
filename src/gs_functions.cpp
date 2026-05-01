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

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "gs_functions.h"

using namespace std;

void sedMATRIX(int* (&step), vector<int>& vec, vector<int> const& pos, int const num, int const N){
  for(auto p: pos){
    vec[p] = num;
  }
  for(int n = 0; n<N; n++){
    step[n*N+num-1] = vec[n];
  }
}

void sedVECTOR(vector<int>& vec, vector<int> const& pos, int const num){
  for(auto p: pos){
    vec[p] = num;
  }
}

double simI(double* const (&W), vector<int> const& res, int const num){
  int N = static_cast<int>(res.size());

  // Collect cluster-member indices in one O(N) pass instead of branching on
  // res[]==num inside an O(N^2) double loop. The visit order over (ro, co)
  // pairs is identical to the original (mem[] is ascending), so the
  // floating-point accumulation order — and thus the bit-exact result — is
  // preserved.
  vector<int> mem;
  mem.reserve(N);
  for(int i = 0; i < N; i++){
    if(res[i] == num) mem.push_back(i);
  }

  int s = static_cast<int>(mem.size());
  if(s <= 1) return 2.0;

  double D = 0;
  for(int i = 0; i < s; i++){
    int ro = mem[i];
    for(int j = 0; j < i; j++){
      int co = mem[j];
      D += 2 * W[ro*N + co];
    }
  }

  return D / (static_cast<double>(s) * (s - 1));
}

int whichMIN(vector<double> const& vec){
  int N      = vec.size();
  int p      = 0;
  // Use +infinity as sentinel so any finite value wins, regardless of how
  // simI is tuned in the future (formerly hard-coded to 100).
  double old = std::numeric_limits<double>::infinity();
  for(int n = 0; n < N; n++){
    // Bug fix: update 'old' only when a new minimum is found,
    // so it tracks the running minimum (not just the previous element)
    if(vec[n] < old){
      p   = n;
      old = vec[n];
    }
  }
  return(p+1);
}

double gev(double const& x, double const& mu){
  double theta = mu*(1-mu)/3.5;
  double gamma = exp(-3*mu)-1;

  if(gamma == 0){
    return mu - theta*log(-log(x));
  }
  else{
    return mu +( pow(-log(x),-gamma)-1 )*theta/gamma;
  }
}
