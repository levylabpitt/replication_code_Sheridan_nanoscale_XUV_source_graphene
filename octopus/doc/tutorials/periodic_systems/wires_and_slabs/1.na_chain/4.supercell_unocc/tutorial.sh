#!/usr/bin/env bash

# these variables need to be defined:
# HELPER_DIR=~/HUGO/octopus-documentation/scripts/
# OCTOPUS_TOP=~/Octopus/octopus/


octopus > log

cp static/bandstructure bandstructure_3D

cp tutorial.sh $OCTOPUS_TOP/doc/tutorials/periodic_systems/wires_and_slabs/1.na_chain/4.supercell_unocc/
