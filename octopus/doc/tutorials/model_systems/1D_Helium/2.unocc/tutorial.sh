#!/usr/bin/env bash

# these variables need to be defined:
# HELPER_DIR=~/HUGO/octopus-documentation/scripts/
# OCTOPUS_TOP=~/Octopus/octopus/

octopus > log

$HELPER_DIR/extract_generic.sh static/info 'Theory Level' 'END   ===' > info.txt

gnuplot plot.gp

cp inp tutorial.sh plot.gp *.eps static/eigenvalues  $OCTOPUS_TOP/doc/tutorials/model_systems/1D_Helium/2.unocc/
