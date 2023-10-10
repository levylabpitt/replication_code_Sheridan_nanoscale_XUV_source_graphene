#!/usr/bin/env bash

# these variables need to be defined:
# HELPER_DIR=~/HUGO/octopus-documentation/scripts/
# OCTOPUS_TOP=~/Octopus/octopus/

cp $OCTOPUS_TOP/testsuite/tutorials/06-octopus_basics-periodic_systems.01-silicon.inp inp

octopus > log

$HELPER_DIR/extract_iter.sh log # creates header.txt, first_iter.txt, last_iter.txt and footer.txt
$HELPER_DIR/extract.sh log 'Space' > Space.txt
$HELPER_DIR/extract.sh log 'Grid' > Grid.txt
$HELPER_DIR/extract.sh log 'Symmetries' > Symmetries.txt
$HELPER_DIR/extract.sh log 'Lattice' > Lattice.txt

$HELPER_DIR/extract_generic.sh log 'Checking if' 'Input:' | head -n -1  > k-points.txt

cp tutorial.sh *.txt *.log log $OCTOPUS_TOP/doc/tutorials/periodic_systems/getting_started_with_periodic_systems/1.start/
