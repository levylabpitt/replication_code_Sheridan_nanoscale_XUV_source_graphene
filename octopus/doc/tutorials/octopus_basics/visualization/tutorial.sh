#!/usr/bin/env bash

# these variables need to be defined:
# HELPER_DIR=~/HUGO/octopus-documentation/scripts/
# OCTOPUS_TOP=~/Octopus/octopus/


octopus > log

gnuplot plot.gp

cp *.sh *.eps $OCTOPUS_TOP/doc/tutorials/octopus_basics/visualization/
