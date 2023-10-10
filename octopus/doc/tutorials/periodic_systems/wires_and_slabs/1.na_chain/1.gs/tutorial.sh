#!/usr/bin/env bash

# these variables need to be defined:
# HELPER_DIR=~/HUGO/octopus-documentation/scripts/
# OCTOPUS_TOP=~/Octopus/octopus/

octopus > log

$HELPER_DIR/extract_iter.sh log # creates header.txt, first_iter.txt, last_iter.txt and footer.txt
$HELPER_DIR/extract.sh log 'Space' > Space.txt
$HELPER_DIR/extract.sh log 'Grid' > Grid.txt
$HELPER_DIR/extract.sh log 'Symmetries' > Symmetries.txt
$HELPER_DIR/extract.sh log 'Lattice' > Lattice.txt
$HELPER_DIR/extract.sh log 'Hartree' > Hartree.txt


cp tutorial.sh *.txt *.log log $OCTOPUS_TOP/doc/tutorials/periodic_systems/wires_and_slabs/1.na_chain/1.gs/
