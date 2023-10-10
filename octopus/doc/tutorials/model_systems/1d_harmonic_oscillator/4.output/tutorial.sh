#!/usr/bin/env bash

# these variables need to be defined:
# HELPER_DIR=~/HUGO/octopus-documentation/scripts/
# OCTOPUS_TOP=~/Octopus/octopus/

octopus > log

gnuplot plot.gp


cp tutorial.sh plot.gp *.eps $OCTOPUS_TOP/doc/tutorials/model_systems/1d_harmonic_oscillator/4.output/
