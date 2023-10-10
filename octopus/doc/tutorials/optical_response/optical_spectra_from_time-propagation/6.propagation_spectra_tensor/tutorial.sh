#!/usr/bin/env bash

# these variables need to be defined:
# HELPER_DIR=~/HUGO/octopus-documentation/scripts/
# OCTOPUS_TOP=~/Octopus/octopus/

rm *.txt *.log

head -n 17 cross_section_tensor > cross_section_tensor-head.txt
gnuplot plot.gp

cp tutorial.sh cross_section_tensor-head.txt *.eps $OCTOPUS_TOP/doc/tutorials/optical_response/optical_spectra_from_time-propagation/6.propagation_spectra_tensor/
