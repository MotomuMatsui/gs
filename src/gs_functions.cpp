/******************************************\
| Graph Splitting Method v2.2 (2018/11/07) |
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
#include <vector>

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

//simI function
double simI(double* const (&W), vector<int> const& res, int const num){
  int N    = res.size();
  int s    = count(res.begin(), res.end(), num);
  double D = 0;

  if(s<=1){
    D = 2;
  }
  else{
    for(int ro=0; ro < N; ro++){
      if(res[ro] == num){
        for(int co=0; co < ro; co++){
          if(res[co] == num){
            D += 2*W[ro*N+co];
          }
        }
      }
    }

    D /= (s*(s-1));
  }

  return D;
}

//Minimum factor
int whichMIN(vector<double> const& vec){
  int N = vec.size();
  int p = 0;
  double old = 100;
  for(int n = 0; n < N; n++){
    p   = (vec[n]<old)? n: p;
    old = vec[n];
  }

  return(p+1);
}

//Generalized Extreme Value function (inverse function)                                                                        
double gev(double const& x, double const& mu){
  double theta = mu*(1-mu)/3;
  double gamma = exp(-3*mu)-1;

  if(gamma == 0){ //Gummbel distribution                                                                                        
    return mu - theta*log(-log(x));
  }
  else{
    return mu +( pow(-log(x),-gamma)-1 )*theta/gamma;
  }
}
