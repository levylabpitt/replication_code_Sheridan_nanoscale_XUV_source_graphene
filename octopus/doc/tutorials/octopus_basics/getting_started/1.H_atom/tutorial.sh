#!/usr/bin/env bash

# these variables need to be defined:
# HELPER_DIR=~/HUGO/octopus-documentation/scripts/
# OCTOPUS_TOP=~/Octopus/octopus/

rm *.log *.txt
octopus > log

$HELPER_DIR/extract_iter.sh log # creates header.txt, first_iter.txt, last_iter.txt and footer.txt
$HELPER_DIR/extract.sh log 'Calculation Mode' > Calculation_mode.txt
$HELPER_DIR/extract.sh log 'Space' > Space.txt
$HELPER_DIR/extract.sh log 'Grid' > Grid.txt
$HELPER_DIR/extract.sh log 'Species' > Species.txt

$HELPER_DIR/extract_generic.sh log 'initial LCAO'  'restart'  > lcao.txt
$HELPER_DIR/extract_generic.sh static/info 'Eigenvalues'  'Dipole' | head -n -1 > info.txt

cp *.txt $OCTOPUS_TOP/doc/tutorials/octopus_basics/getting_started/1.H_atom/
