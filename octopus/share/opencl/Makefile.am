## Process this file with automake to produce Makefile.in

## Copyright (C) 2010 X. Andrade
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
## 02110-1301, USA.
##

opencldir = $(pkgdatadir)/opencl

dist_opencl_DATA =              \
	cl_global.h             \
	cl_complex.h            \
	cl_rtype.h              \
	cl_reduce.h             \
	cuda_compat.h           \
	vpsi.cl                 \
	set_zero.cl             \
	set_one.cl              \
	operate.cl              \
	axpy.cl                 \
	projector.cl            \
	pack.cl                 \
	copy.cl                 \
	mesh_batch_single.cl	\
	density.cl              \
	phase.cl                \
	phase_spiral.cl         \
	mesh_batch.cl           \
	rotate.cl               \
	points.cl               \
	trsm.cl                 \
	mul.cl                  \
	mesh_to_cube.cl         \
	convert.cl              \
	boundaries.cl           \
	ghost.cl                \
	forces.cl               \
	mod_sqr.cl              \
	uvw_to_xyz.cl           \
	curl.cl                 \
	pml.cl                  \
	dftu_projector.cl

# In order to get `make check' working with a `VPATH' build (e. g. to use
# the `distcheck' target) the pseudopotential files have to be accessible
# in the build tree.
check:
	for f in $(dist_opencl_DATA); do \
	  if [ ! -e $$f -o $$f -ot $(srcdir)/$$f ]; then \
	    cp $(srcdir)/$$f $$f; \
	  fi \
	done

check-short: check
check-long: check

clean-local:
	if [ -n "$(top_builddir)" -a "$(top_builddir)" != "$(top_srcdir)" ]; then \
	  rm -f $(dist_opencl_DATA); \
	fi

.PHONY: check check-short check-long
