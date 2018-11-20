/******************************************\
| Graph Splitting Method v2.3 (2018/11/16) |
|                                          |
| Copyright (c) 2015-2018 Motomu Matsui    |
|     Distributed under the GNU GPL        |
|                                          |
|     Matsui M and Iwasaki W (2018)        |
|     Systematic Biology, xx:xx-xx.        |
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
