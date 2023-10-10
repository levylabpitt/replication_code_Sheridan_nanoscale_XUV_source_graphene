#!/usr/bin/env bash

# these variables need to be defined:
# HELPER_DIR=~/HUGO/octopus-documentation/scripts/
# OCTOPUS_TOP=~/Octopus/octopus/

rm *.txt
oct-center-geom > log

$HELPER_DIR/extract_generic.sh log ' Symmetries '  'Parser warning' > Symmetries.txt

cp tutorial.sh adjusted.xyz *.txt $OCTOPUS_TOP/doc/tutorials/octopus_basics/centering_a_geometry/
