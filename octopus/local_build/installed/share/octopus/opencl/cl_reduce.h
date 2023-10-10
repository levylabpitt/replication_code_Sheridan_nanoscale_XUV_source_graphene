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

// Define the CUDA warpReduce function

#ifdef CUDA
// shuffle instructions wrappers
#if __CUDACC_VER_MAJOR__ >= 9

#define MASK_ALL_WARP 0xFFFFFFFF

#define warpShflDown(var, delta)   __shfl_down_sync (MASK_ALL_WARP, var, delta)

#else

#define warpShflDown(var, delta)   __shfl_down (var, delta)

#endif

// See for instance 
// https://developer.nvidia.com/blog/faster-parallel-reductions-kepler/

__device__ inline double dwarpReduce(double val)
{
#pragma unroll
  for (int offset = warpSize/2; offset > 0; offset /= 2){
    val += warpShflDown(val, offset);
  }
  return val;
}

__device__ inline double2 zwarpReduce(double2 val)
{
#pragma unroll
  for (int offset = warpSize/2; offset > 0; offset /= 2){
    val.x += warpShflDown(val.x, offset);
    val.y += warpShflDown(val.y, offset);
  }
  return val;
}

#endif



/*
 Local Variables:
 mode: c
 coding: utf-8
 End:
*/
