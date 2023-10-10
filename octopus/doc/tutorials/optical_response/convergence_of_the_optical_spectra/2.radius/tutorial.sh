#!/usr/bin/env bash

# these variables need to be defined:
# HELPER_DIR=~/HUGO/octopus-documentation/scripts/
# OCTOPUS_TOP=~/Octopus/octopus/

. radius.sh

gnuplot plot.gp

cp cross_section_vector-6.5  $OCTOPUS_TOP/doc/tutorials/optical_response/optical_spectra_from_casida/4.spectrum/
cp cross_section_vector-6.5  $OCTOPUS_TOP/doc/tutorials/optical_response/triplet_excitations/1.time_propagation/3.spectrum/cross_section_vector-singlets

cp radius.sh tutorial.sh  Absorption_spectrum_CH4_radius.eps $OCTOPUS_TOP/doc/tutorials/optical_response/convergence_of_the_optical_spectra/2.radius/
