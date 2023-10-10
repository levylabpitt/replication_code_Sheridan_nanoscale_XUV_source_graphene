#!/usr/bin/env bash

# these variables need to be defined:
# HELPER_DIR=~/HUGO/octopus-documentation/scripts/
# OCTOPUS_TOP=~/Octopus/octopus/

octopus > log

$HELPER_DIR/extract.sh log 'Hartree' > Hartree.txt


cp tutorial.sh Hartree.txt $OCTOPUS_TOP/doc/tutorials/periodic_systems/wires_and_slabs/2.h-bn_monolayer/1.gs/
