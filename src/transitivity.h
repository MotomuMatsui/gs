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

#ifndef TRANSITIVITY_H
#define TRANSITIVITY_H

#include <cblas.h>

double transitivity(double* const (&oW), int const& size);

#endif /*TRANSITIVITY_H*/
