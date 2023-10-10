#!/usr/bin/env bash

# these variables need to be defined:
# HELPER_DIR=~/HUGO/octopus-documentation/scripts/
# OCTOPUS_TOP=~/Octopus/octopus/


oct-dielectric-function > dielectric.log

gnuplot plot.gp


cp tutorial.sh plot.gp *.eps $OCTOPUS_TOP/doc/tutorials/periodic_systems/optical_spectra_of_solids_gauge_field/3.spectra
