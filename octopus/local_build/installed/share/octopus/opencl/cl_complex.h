/*
 Copyright (C) 2012 X. Andrade

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2, or (at your option)
 any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 02110-1301, USA.

*/

#ifndef __CL_COMPLEX_H__
#define __CL_COMPLEX_H__

#include <cl_global.h>

inline double2 complex_number(const double x, const double y){
#ifdef CUDA
  return double2(x, y);
#else
    return (double2) (x, y);
#endif
}

inline double2 complex_conj(const double2 a){
#ifdef CUDA
  return double2(a.x, -a.y);
#else
    return (double2) (a.x, -a.y);
#endif
}

inline double2 complex_mul(const double2 a, const double2 b){
#ifdef CUDA
  return double2(a.x*b.x - a.y*b.y, a.y*b.x + a.x*b.y);
#else
  return (double2)(a.x*b.x - a.y*b.y, a.y*b.x + a.x*b.y);
#endif
}

inline double2 complex_dotp(const double2 a, const double2 b){
#ifdef CUDA
  return double2(a.x*b.x + a.y*b.y,-a.y*b.x + a.x*b.y);
#else
  return (double2)(a.x*b.x + a.y*b.y,-a.y*b.x + a.x*b.y);
#endif
}


inline double2 complex_div(const double2 a, const double2 b){
  double2 c = b*b;
  return complex_mul(a, complex_conj(b))/(c.x + c.y);
}

inline double complex_real(const double2 a){
  return a.x;
}

inline double complex_imag(const double2 a){
  return a.y;
}


#endif

/*
 Local Variables:
 mode: c
 coding: utf-8
 End:
*/
