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

#ifndef GS_FUNCTION_H
#define GS_FUNCTION_H

#include <algorithm>
#include <cmath>
#include <vector>

void sedMATRIX(int* (&step), std::vector<int>& vec, std::vector<int> const& pos, int const num, int const N);
void sedVECTOR(std::vector<int>& vec, std::vector<int> const& pos, int const num);
double simI(double* const (&W), std::vector<int> const& res, int const num);
int whichMIN(std::vector<double> const& vec);
double gev(double const& x, double const& mu);

#endif /*GS_FUNCTION_H*/
