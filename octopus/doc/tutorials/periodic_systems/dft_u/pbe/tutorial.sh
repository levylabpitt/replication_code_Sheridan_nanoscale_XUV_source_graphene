#!/usr/bin/env bash

# these variables need to be defined:
# HELPER_DIR=~/HUGO/octopus-documentation/scripts/
# OCTOPUS_TOP=~/Octopus/octopus/


octopus > log

$HELPER_DIR/extract_iter.sh log # creates header.txt, first_iter.txt, last_iter.txt and footer.txt
$HELPER_DIR/extract.sh log 'Symmetries' > Symmetries.txt
$HELPER_DIR/extract.sh log 'Lattice' > Lattice.txt

tail -n 40 log > end 
$HELPER_DIR/extract_generic.sh end 'Total Magnetic Moment' 'Elapsed' | head -n -2  > Moments.txt

grep -A2 'Direct' static/info > Gaps.txt
$HELPER_DIR/extract_generic.sh log 'Checking if' 'Input:' | head -n -1  > k-points.txt


cp tutorial.sh *.txt *.log log $OCTOPUS_TOP/doc/tutorials/periodic_systems/dft_u/pbe/
