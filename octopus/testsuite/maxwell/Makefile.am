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

sharedir = $(pkgdatadir)/testsuite/maxwell

dist_share_DATA =                                                                   \
	01-free-propagation.test                                                        \
	01-free-propagation.01-1_pulse_td.inp                                           \
    01-free-propagation.02-2_pulses_td.inp                                          \
    02-external-current.test                                                        \
    02-external-current.01-gaussian_current_pulse.inp                               \
    02-external-current.02-gaussian_current_pulse_with_pml.inp                      \
    02-external-current.03-gaussian_current_pulse_with_mask.inp                     \
    03-linear-medium.test                                                           \
	03-linear-medium.01-cosinoidal_pulse_td.inp                                     \
	03-linear-medium.02-cosinoidal_pulse_td_pml.inp                                 \
	03-linear-medium.03-cosinoidal_pulse_td_pml_unpacked.inp                        \
    04-linear-medium-from-file.test                                                 \
	04-linear-medium-from-file.01-cube_medium_from_file.inp                         \
	04-linear-medium-from-file.01-cube.off                                          \
    05-plane_waves.test                                                             \
    05-plane_waves.01-pulse_pml.inp                                                 \
    06-circular-polarization.test                                                   \
    06-circular-polarization.01-1_puse_circ_pol.inp                                 \
	07-mode-resolved-maxwell-ks-propagation.test                                    \
	07-mode-resolved-maxwell-ks-propagation.01-1D-Helium-ground-state.inp           \
	07-mode-resolved-maxwell-ks-propagation.02-kick-of-the-electronic-subsystem.inp \
	08-restart-maxwell.test                                                         \
	08-restart-maxwell.01-free-propagation_fromscratch.inp                          \
	08-restart-maxwell.01-free-propagation_restarted.inp                            \
	08-restart-maxwell.02-external-current-pml_fromscratch.inp                      \
	08-restart-maxwell.02-external-current-pml_restarted.inp                        \
    08-restart-maxwell.03-linear-medium_fromscratch.inp                             \
	08-restart-maxwell.03-linear-medium_restarted.inp                               \
    08-restart-maxwell.03-linear-medium-cube.off                                    \
    09-drude-medium-from-file.test                                                  \
    09-drude-medium-from-file.01-gold-np.inp                                        \
	09-drude-medium-from-file.01-gold-np-r80nm.off                                  \
    10-current-to-maxwell.test                                                      \
    10-current-to-maxwell.01-all-forces-non-self-consistent.inp

CLEANFILES = *~ *.bak
