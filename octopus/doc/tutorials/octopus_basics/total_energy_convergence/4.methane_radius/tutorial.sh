#!/usr/bin/env bash

# these variables need to be defined:
# HELPER_DIR=~/HUGO/octopus-documentation/scripts/
# OCTOPUS_TOP=~/Octopus/octopus/

rm radius.log
. ./radius.sh

gnuplot plot.gp

cp *.sh *.log *.eps $OCTOPUS_TOP/doc/tutorials/octopus_basics/total_energy_convergence/4.methane_radius/
