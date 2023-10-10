#!/usr/bin/env bash

# these variables need to be defined:
# HELPER_DIR=~/HUGO/octopus-documentation/scripts/
# OCTOPUS_TOP=~/Octopus/octopus/

rm *.log *.txt
export OMP_NUM_THREADS=2
mpirun -n 1 octopus > log

$HELPER_DIR/extract_iter.sh log # creates header.txt, first_iter.txt, last_iter.txt and footer.txt
$HELPER_DIR/extract.sh log 'Parallelization' > parallelization.txt
$HELPER_DIR/extract_generic.sh log 'Info: Parallelization' 'Info: Generating' | head -n -1 > parstates.txt


cp inp tutorial.sh parallelization.txt parstates.txt $OCTOPUS_TOP/doc/tutorials/hpc/parallelization/4.omp/
