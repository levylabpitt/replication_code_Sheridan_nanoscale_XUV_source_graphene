#!/usr/bin/env bash

# these variables need to be defined:
# HELPER_DIR=~/HUGO/octopus-documentation/scripts/
# OCTOPUS_TOP=~/Octopus/octopus/

rm *.log *.txt
mpirun -n 2 octopus > log

$HELPER_DIR/extract_iter.sh log # creates header.txt, first_iter.txt, last_iter.txt and footer.txt
$HELPER_DIR/extract.sh log 'Parallelization' > parallelization.txt
$HELPER_DIR/extract_generic.sh log 'Info: Parallelization' 'Info: Generating' | head -n -1 > parstates.txt


cp tutorial.sh parallelization.txt parstates.txt $OCTOPUS_TOP/doc/tutorials/hpc/parallelization/2.td/
