#!/usr/bin/env bash

# these variables need to be defined:
# HELPER_DIR=~/HUGO/octopus-documentation/scripts/
# OCTOPUS_TOP=~/Octopus/octopus/

rm *.txt
octopus > log

$HELPER_DIR/extract.sh log 'Symmetries' > Symmetries.txt


cp tutorial.sh *.txt $OCTOPUS_TOP/doc/tutorials/periodic_systems/optical_spectra_of_solids_conductivity/1.gs
