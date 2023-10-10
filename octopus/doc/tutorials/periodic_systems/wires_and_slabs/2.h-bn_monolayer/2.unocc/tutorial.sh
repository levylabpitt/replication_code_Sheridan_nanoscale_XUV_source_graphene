#!/usr/bin/env bash

# these variables need to be defined:
# HELPER_DIR=~/HUGO/octopus-documentation/scripts/
# OCTOPUS_TOP=~/Octopus/octopus/

octopus > log

gnuplot plot.gp


cp tutorial.sh plot.gp *.eps $OCTOPUS_TOP/doc/tutorials/periodic_systems/wires_and_slabs/2.h-bn_monolayer/2.unocc/
