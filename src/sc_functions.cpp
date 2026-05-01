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
#include <functional>
#include <iterator>
#include <map>
#include <memory>
#include <utility>
#include <vector>

#include "sc_functions.h"

using namespace std;

void subMATRIX(double* const (&Wo), double* (&W), double* (&D), vector<int> const& res, int const num){

  int N = static_cast<int>(res.size());

  // Collect cluster member positions in one O(N) pass, then loop only over
  // members.  This eliminates two layers of `res[i]==num` branching from the
  // hot inner loop.  D's accumulation order over (rn,cn) pairs is identical
  // to the original (mem[] is monotonically ascending), so the floating-point
  // result is bit-exact.
  vector<int> mem;
  mem.reserve(N);
  for(int i = 0; i < N; i++){
    if(res[i] == num) mem.push_back(i);
  }
  int s = static_cast<int>(mem.size());

  for(int rn = 0; rn < s; rn++){
    int ro = mem[rn];
    const double* Wrow = Wo + static_cast<size_t>(ro) * N;
    for(int cn = 0; cn < rn; cn++){
      int co     = mem[cn];
      double p   = Wrow[co];
      W[rn*s+cn] = p;
      W[cn*s+rn] = p;
      D[rn]     += p;
      D[cn]     += p;
    }
    double p   = Wrow[ro];
    W[rn*s+rn] = p;
    D[rn]     += p;
  }

  for(int p = 0; p<s; p++){
    D[p] = 1/sqrt(D[p]);
  }
}

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

//Best cut position
typedef pair<double, int> P;
// Internal helper — `static` gives internal linkage to prevent ODR clashes
// across translation units.
static bool comp(const P &a, const P &b){
  return a.first > b.first;
}

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

void whichCUT(double* const (&z), int const col, double* const (&D), int*& qi, int& cut, int const N){
  vector<P> pairs(N);
  auto qs = unique_ptr<double[]>(new double[N]());

  //Inner product: <eigenvector,D>
  for(int p = 0; p < N; p++){
    pairs[p] = make_pair(z[col+p*N]*D[p], p);
  }

  //Sorting values
  sort(pairs.begin(), pairs.end(), comp);
  for(int p = 0; p < N; p++){
    qs[p]  = pairs[p].first;
    qi[p]  = pairs[p].second;
  }

  cut = 0;
  double best = 0;
  for(int n = 0; n < N-1; n++){
    auto def = abs(qs[n] - qs[n+1]);
    if(def > best){
      cut  = n;
      best = def;
    }
  }
}

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

void splitVECTOR(vector<int> const& res, int const num, int* const (&qi), int const cut, vector<int>& a, vector<int>& b, int const N){

  vector<int> index;
  int n = 0;
  for(int v: res){
    if(v == num){
      index.push_back(n);
    }
    n ++;
  }

  for(int p = 0; p <=cut; p++){
    a.push_back(index[qi[p]]);
  }
  for(int p = cut+1; p < N; p++){
    b.push_back(index[qi[p]]);
  }
}

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
