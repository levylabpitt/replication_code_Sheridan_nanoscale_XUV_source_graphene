#!/usr/bin/env bash

# these variables need to be defined:
# HELPER_DIR=~/HUGO/octopus-documentation/scripts/
# OCTOPUS_TOP=~/Octopus/octopus/


gnuplot plot.gp

cp tutorial.sh plot.gp *.eps $OCTOPUS_TOP/doc/tutorials/periodic_systems/wires_and_slabs/1.na_chain/5.plot/
