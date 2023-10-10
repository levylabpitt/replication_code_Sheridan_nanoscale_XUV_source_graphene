#!/usr/bin/env bash

# these variables need to be defined:
# HELPER_DIR=~/HUGO/octopus-documentation/scripts/
# OCTOPUS_TOP=~/Octopus/octopus/

octopus > log

$HELPER_DIR/extract_iter.sh log 
$HELPER_DIR/extract.sh log 'Species' > species.txt
$HELPER_DIR/extract.sh log 'Theory Level' > theory_level.txt


cp tutorial.sh species.txt theory_level.txt last_iter.txt  $OCTOPUS_TOP/doc/tutorials/model_systems/1d_harmonic_oscillator/1.main/
