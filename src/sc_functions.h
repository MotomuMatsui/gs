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

#ifndef SC_FUNCTIONS_H
#define SC_FUNCTIONS_H

#include <algorithm>
#include <cmath>
#include <functional>
#include <iterator>
#include <map>
#include <utility>
#include <vector>

typedef std::pair<double, int> P;
// Note: the comparator used by whichCUT lives in sc_functions.cpp with internal
// linkage (static) to avoid ODR clashes; intentionally not declared here.

void subMATRIX(double* const (&Wo), double* (&W), double* (&D), std::vector<int> const& res, int const num);
void whichCUT(double* const (&z), double* const (&D), int*& qi, int& cut, int const N);
void splitVECTOR(std::vector<int> const& res, int const num, int* const (&qi), int const cut, std::vector<int>& a, std::vector<int>& b, int const N);

#endif /*SC_FUNCTIONS_H*/
