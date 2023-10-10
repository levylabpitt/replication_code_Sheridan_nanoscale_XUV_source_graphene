#!/usr/bin/env bash

# these variables need to be defined:
# HELPER_DIR=~/HUGO/octopus-documentation/scripts/
# OCTOPUS_TOP=~/Octopus/octopus/

. ./spacing.sh 

gnuplot plot.gp

cp tutorial.sh  *.eps $OCTOPUS_TOP/doc/tutorials/optical_response/convergence_of_the_optical_spectra/1.spacing/
