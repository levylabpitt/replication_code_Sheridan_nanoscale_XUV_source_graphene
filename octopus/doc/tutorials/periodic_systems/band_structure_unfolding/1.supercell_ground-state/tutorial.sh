#!/usr/bin/env bash

# these variables need to be defined:
# HELPER_DIR=~/HUGO/octopus-documentation/scripts/
# OCTOPUS_TOP=~/Octopus/octopus/

octopus > log

cp tutorial.sh $OCTOPUS_TOP/doc/tutorials/periodic_systems/band_structure_unfolding/1.supercell_ground-state/
