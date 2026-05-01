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

#ifndef GS_H
#define GS_H

#include <tuple>
#include <vector>

#include "sc.h"
#include "gs_functions.h"

void GS(double* const (&W), int* (&step), int const& size);

#endif /*GS_H*/
