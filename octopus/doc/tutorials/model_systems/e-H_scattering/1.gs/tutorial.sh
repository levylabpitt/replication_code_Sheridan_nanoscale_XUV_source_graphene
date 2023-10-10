#!/usr/bin/env bash

# these variables need to be defined:
# HELPER_DIR=~/HUGO/octopus-documentation/scripts/
# OCTOPUS_TOP=~/Octopus/octopus/

octopus > log

$HELPER_DIR/extract.sh log "Theory Level" > theory_level.txt

cp tutorial.sh theory_level.txt $OCTOPUS_TOP/doc/tutorials/model_systems/e-H_scattering/1.gs/
