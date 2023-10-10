#!/usr/bin/env bash

# these variables need to be defined:
# HELPER_DIR=~/HUGO/octopus-documentation/scripts/
# OCTOPUS_TOP=~/Octopus/octopus/

oct-unfold > unfold-run.log

head -n 15 static/ake.dat > ake.txt

gnuplot plot.gp

cp tutorial.sh ake.txt plot.gp  *.eps $OCTOPUS_TOP/doc/tutorials/periodic_systems/band_structure_unfolding/4.unfolding/
