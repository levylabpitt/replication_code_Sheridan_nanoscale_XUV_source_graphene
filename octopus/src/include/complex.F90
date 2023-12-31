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

#define R_TCOMPLEX 1

#define R_TYPE      CMPLX
#define R_BASE      FLOAT
#define R_DOUBLE    complex(r8)
#define R_MPITYPE   MPI_CMPLX
#define R_TYPE_VAL  TYPE_CMPLX
#define R_TYPE_CL   'RTYPE_COMPLEX'
#define R_TYPE_IOBINARY TYPE_DOUBLE_COMPLEX
#define R_TOTYPE(x) cmplx(x, M_ZERO, REAL_PRECISION)

#define R_CONJ(x)   conjg(x)
#define R_REAL(x)   real(x)
#define R_AIMAG(x)  aimag(x)

#define R_SIZEOF    16
#define R_ADD       2
#define R_MUL       6

#ifdef __GFORTRAN__
#define X(x)        z/**/x
#define pX(x)       pz/**/x
#define aX(x,y)     x/**/z/**/y
#else
#define X(x)        z ## x
#define pX(x)       pz ## x
#define aX(x,y)     x ## z ## y
#endif


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
