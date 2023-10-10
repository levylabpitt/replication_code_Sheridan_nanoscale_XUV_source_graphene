#!/usr/bin/env bash

# these variables need to be defined:
# HELPER_DIR=~/HUGO/octopus-documentation/scripts/
# OCTOPUS_TOP=~/Octopus/octopus/

octopus > log

$HELPER_DIR/extract_iter.sh log 


cp tutorial.sh last_iter.txt $OCTOPUS_TOP/doc/tutorials/model_systems/1d_harmonic_oscillator/2.lcao_start/
