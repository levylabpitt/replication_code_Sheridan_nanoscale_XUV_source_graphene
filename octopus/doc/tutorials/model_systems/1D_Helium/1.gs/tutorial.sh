#!/usr/bin/env bash

# these variables need to be defined:
# HELPER_DIR=~/HUGO/octopus-documentation/scripts/
# OCTOPUS_TOP=~/Octopus/octopus/

rm -rf restart/
octopus > log

$HELPER_DIR/extract_generic.sh static/info 'Theory Level' 'END   ===' > info.txt

cp tutorial.sh info.txt  $OCTOPUS_TOP/doc/tutorials/model_systems/1D_Helium/1.gs/
