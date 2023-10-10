#!/usr/bin/env bash

# these variables need to be defined:
# HELPER_DIR=~/HUGO/octopus-documentation/scripts/
# OCTOPUS_TOP=~/Octopus/octopus/

oct-casida_spectrum > spectrum.log 

cp casida/spectrum.casida spectrum.casida-triplets
cp $OCTOPUS_TOP/doc/tutorials/optical_response/optical_spectra_from_casida/4.spectrum/spectrum.casida-singlets .
cp $OCTOPUS_TOP/doc/tutorials/optical_response/triplet_excitations/1.time_propagation/3.spectrum/cross_section_vector-triplets cross_section_vector
gnuplot plot.gp

cp tutorial.sh plot.gp Triplet_casida_td_CH4.eps $OCTOPUS_TOP/doc/tutorials/optical_response/triplet_excitations/3.comparison/
