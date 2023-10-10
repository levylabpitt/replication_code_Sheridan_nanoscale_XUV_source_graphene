#!/usr/bin/env bash

# these variables need to be defined:
# HELPER_DIR=~/HUGO/octopus-documentation/scripts/
# OCTOPUS_TOP=~/Octopus/octopus/

octopus > log
head -5 static/bandstructure > bandstructure-head.txt

gnuplot plot.gp

cp tutorial.sh plot.gp *.eps *.txt $OCTOPUS_TOP/doc/tutorials/periodic_systems/getting_started_with_periodic_systems/4.band-structure
