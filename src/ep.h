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

#ifndef EP_H
#define EP_H

#include <algorithm>
#include <cmath> 
#include <functional>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "sc.h"
#include "gs_functions.h"
#include "format.h"

using namespace std;

void EP_fbs(double* const (&W), unordered_map<string, double>& ep, function<double()>& R, int const& size);
void EP_tbe(double* const (&W), int* const (&list_ori), unordered_map<string, double>& ep, function<double()>& R, int const& size);

#endif /*EP_H*/
