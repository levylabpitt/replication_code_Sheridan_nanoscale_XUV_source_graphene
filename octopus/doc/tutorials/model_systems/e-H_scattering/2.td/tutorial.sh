#!/usr/bin/env bash

# these variables need to be defined:
# HELPER_DIR=~/HUGO/octopus-documentation/scripts/
# OCTOPUS_TOP=~/Octopus/octopus/

octopus > log

gnuplot plot.gp

./convert.sh

cp tutorial.sh theory_level.txt $OCTOPUS_TOP/doc/tutorials/model_systems/e-H_scattering/2.td/
cp scattering.gif $OCTOPUS_TOP/doc/html/images/
