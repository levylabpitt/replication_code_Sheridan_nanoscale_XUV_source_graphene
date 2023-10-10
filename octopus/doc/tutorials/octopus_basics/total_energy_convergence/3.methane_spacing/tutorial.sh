#!/usr/bin/env bash

# these variables need to be defined:
# HELPER_DIR=~/HUGO/octopus-documentation/scripts/
# OCTOPUS_TOP=~/Octopus/octopus/

rm spacing.log
. ./spacing.sh

gnuplot plot.gp

cp *.sh *.log *.eps *.txt $OCTOPUS_TOP/doc/tutorials/octopus_basics/total_energy_convergence/3.methane_spacing/
