#!/usr/bin/env bash

# these variables need to be defined:
# HELPER_DIR=~/HUGO/octopus-documentation/scripts/
# OCTOPUS_TOP=~/Octopus/octopus/

oct-unfold > unfold.log

cp tutorial.sh unfold*.dat $OCTOPUS_TOP/doc/tutorials/periodic_systems/band_structure_unfolding/2.unfolding_setup/
