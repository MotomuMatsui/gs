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

#ifndef SC_FUNCTIONS_H
#define SC_FUNCTIONS_H

#include <algorithm>
#include <cmath> 
#include <functional>
#include <iterator>
#include <map>
#include <utility>
#include <vector>

using namespace std;

typedef pair<double, int> P;
bool comp(const P &a, const P &b);
void subMATRIX(double* const (&Wo), double* (&W), double* (&D), vector<int> const& res, int const num);
void whichCUT(double* const (&z), int const col, double* const (&D), int*& qi, int& cut, int const N);
void splitVECTOR(vector<int> const& res, int const num, int* const (&qi), int const cut, vector<int>& a, vector<int>& b, int const N);

#endif /*SC_FUNCTIONS_H*/
