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
#include <algorithm>
#include <cmath> 
#include <functional>
#include <utility>
#include <map>
#include <iterator>

using namespace std;

void subMATRIX(double* const (&Wo), double* (&W), double* (&D), vector<int> const& res, int const num){

  int N = res.size();
  int s = count(res.begin(), res.end(), num);

  int rn   = 0;
  int cn   = 0;
  double p;
  for(int ro=0; ro < N; ro++){
    if(res[ro] == num){
      cn=0;
      for(int co=0; co < ro; co++){
        if(res[co] == num){
          p = Wo[ro*N+co];
          W[rn*s+cn] = p;
          W[cn*s+rn] = p;
          D[rn]     += p;
          D[cn]     += p;
          cn++;
        }
      }
      p = Wo[ro*N+ro];
      W[rn*s+rn] = p;
      D[rn]     += p;

      rn++;
    }
  }

  for(int p = 0; p<s; p++){
    D[p] = 1/sqrt(D[p]);
  }
}

//Best cut position
typedef pair<double, int> P;
bool comp(const P &a, const P &b){
  return a.first > b.first;
}

void whichCUT(double* const (&z), int const col, double* const (&D), int*& qi, int& cut, int const N){
  vector<P> pairs(N);
  double* qs = new double[N]();

  //Inner product: <eigenvector,D>
  for(int p = 0; p < N; p++){
    pairs[p] = make_pair(z[col+p*N]*D[p], p);
  }

  //Sorting values
  sort(pairs.begin(), pairs.end(), comp);
  for(int p = 0; p < N; p++){
    qs[p] = pairs[p].first;
    qi[p] = pairs[p].second;
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

  delete[] qs;
}

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
