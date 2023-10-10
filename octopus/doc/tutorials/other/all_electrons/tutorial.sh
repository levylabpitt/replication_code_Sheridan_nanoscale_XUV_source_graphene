#!/usr/bin/env bash

# these variables need to be defined:
# HELPER_DIR=~/HUGO/octopus-documentation/scripts/
# OCTOPUS_TOP=~/Octopus/octopus/

rm *.log *.txt
octopus > log

grep -A 10 'Eigenvalues \[H\]' static/info > info.txt

. ./spacing.sh

gnuplot plot.gp

cp tutorial.sh *.txt *.eps $OCTOPUS_TOP/doc/tutorials/other/all_electrons/
