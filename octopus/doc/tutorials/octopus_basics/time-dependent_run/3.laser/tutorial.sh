#!/usr/bin/env bash

# these variables need to be defined:
# HELPER_DIR=~/HUGO/octopus-documentation/scripts/
# OCTOPUS_TOP=~/Octopus/octopus/

octopus > log

gnuplot plot.gp

cp tutorial.sh *.eps log $OCTOPUS_TOP/doc/tutorials/octopus_basics/time-dependent_run/3.laser
