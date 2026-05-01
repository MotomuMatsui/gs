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

#ifndef SC_H
#define SC_H

#include <algorithm>
#include <tuple>
#include <vector>

#include "eigen.h"
#include "sc_functions.h"

std::tuple<std::vector<int>,std::vector<int>> spectral_clustering(
    double* const& oW, std::vector<int> const& res, int const& num);

#endif /*SC_H*/
