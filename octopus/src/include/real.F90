!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#define R_TREAL     1

#define R_TYPE      FLOAT
#define R_BASE      FLOAT
#define R_DOUBLE    real(r8)
#define R_MPITYPE   MPI_FLOAT
#define R_TYPE_VAL  TYPE_FLOAT
#define R_TYPE_CL   'RTYPE_DOUBLE'
#define R_TYPE_IOBINARY TYPE_DOUBLE
#define R_TOTYPE(x) real(x, REAL_PRECISION)

#define R_CONJ(x)   (x)
#define R_REAL(x)   (x)
#define R_AIMAG(x)  (M_ZERO)

#ifdef __GFORTRAN__
#define X(x)        d/**/x
#define pX(x)       pd/**/x
#define aX(x,y)     x/**/d/**/y
#else
#define X(x)        d ## x
#define pX(x)       pd ## x
#define aX(x,y)     x ## d ## y
#endif

#define R_SIZEOF    8
#define R_ADD       1
#define R_MUL       1

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
