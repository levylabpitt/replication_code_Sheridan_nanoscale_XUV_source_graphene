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

#ifndef __CL_RTYPE_H__
#define __CL_RTYPE_H__

#include <cl_global.h>

#if defined(RTYPE_DOUBLE)

typedef double rtype;
#ifdef __GFORTRAN__
#define X(x)        d/**/x
#else
#define X(x)        d ## x
#endif
#define MUL(x, y)   ((x)*(y))
#define CONJ(x)     (x)
#define REAL(x)     (x)
#define IMAG(x)     (0.0)

#elif defined(RTYPE_COMPLEX)

typedef double2 rtype;
#ifdef __GFORTRAN__
#define X(x)        z/**/x
#else
#define X(x)        z ## x
#endif
#define MUL(x, y)   complex_mul(x, y)
#define CONJ(x)     complex_conj(x)
#define REAL(x)     complex_real(x)
#define IMAG(x)     complex_imag(x)

#else
#error Type not defined
#endif

#endif

/*
 Local Variables:
 mode: c
 coding: utf-8
 End:
*/
