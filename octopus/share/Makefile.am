## Process this file with automake to produce Makefile.in

## Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

SUBDIRS = pseudopotentials recipes util opencl

dta1dir = $(pkgdatadir)
dist_dta1_DATA = logo                 \
	filter_mask.data              \
	gequad.data                   \
        gdlib.png

dta2dir = $(pkgdatadir)
dist_dta2_DATA = variables varinfo

CLEANFILES = variables varinfo varinfo_orig

varinfo:
	$(top_srcdir)/scripts/mk_varinfo.pl -s $(top_srcdir) -b $(top_builddir)

variables:
	$(top_srcdir)/scripts/mk_varinfo.pl -s $(top_srcdir) -b $(top_builddir)

check:
	@for f in $(dist_dta1_DATA) $(dist_dta2_DATA); do \
	  if [ ! -e $$f -o $$f -ot $(srcdir)/$$f ]; then \
	    cp $(srcdir)/$$f $$f; \
	  fi \
	done

check-short: check
check-long: check

clean-local:
	if [ -n "$(top_builddir)" -a "$(top_builddir)" != "$(top_srcdir)" ]; then \
	  rm -f $(dist_dta1_DATA); \
	fi

.PHONY: check check-short check-long
