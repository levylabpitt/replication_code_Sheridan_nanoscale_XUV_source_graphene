!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!! Copyright (C) 2021 S. Ohlmann
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

#define R_TINTEGER8  1

#ifdef __GFORTRAN__
#define X(x)        l/**/x
#else
#define X(x)        l ## x
#endif

#define R_TYPE      integer(i8)
#define R_BASE      integer(i8)
#define R_TYPE_VAL  TYPE_INTEGER8
#define R_MPITYPE   MPI_INTEGER8
#define R_TYPE_IOBINARY TYPE_INT_64
#define R_TOTYPE(x) (x)

#define R_SIZEOF    8

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
