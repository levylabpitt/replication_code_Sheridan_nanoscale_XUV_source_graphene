#!/usr/bin/env bash

# these variables need to be defined:
# HELPER_DIR=~/HUGO/octopus-documentation/scripts/
# OCTOPUS_TOP=~/Octopus/octopus/

octopus > wftn.log

cp tutorial.sh $OCTOPUS_TOP/doc/tutorials/periodic_systems/band_structure_unfolding/3.obtaining_the_wavefunctions/
