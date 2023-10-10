#!/usr/bin/env bash

# these variables need to be defined:
# HELPER_DIR=~/HUGO/octopus-documentation/scripts/
# OCTOPUS_TOP=~/Octopus/octopus/

rm *.log *.txt
octopus > log

$HELPER_DIR/extract_iter.sh log # creates header.txt, first_iter.txt, last_iter.txt and footer.txt
$HELPER_DIR/extract.sh log 'Parallelization' > parallelization.txt


cp *.txt $OCTOPUS_TOP/doc/tutorials/hpc/parallelization/1.gs/
