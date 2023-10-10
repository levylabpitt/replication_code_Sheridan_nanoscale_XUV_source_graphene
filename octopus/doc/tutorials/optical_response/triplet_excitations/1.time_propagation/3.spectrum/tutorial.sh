#!/usr/bin/env bash

# these variables need to be defined:
# HELPER_DIR=~/HUGO/octopus-documentation/scripts/
# OCTOPUS_TOP=~/Octopus/octopus/

oct-propagation_spectrum > spectrum.log

cp cross_section_vector cross_section_vector-triplets
head -n 27 cross_section_vector > cross_section_vector.txt

gnuplot plot.gp

cp tutorial.sh *.txt plot.gp cross_section_vector-triplets Singlet_triplet_spectrum_CH4.eps $OCTOPUS_TOP/doc/tutorials/optical_response/triplet_excitations/1.time_propagation/3.spectrum/
