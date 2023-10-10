#!/usr/bin/env bash

# these variables need to be defined:
# HELPER_DIR=~/HUGO/octopus-documentation/scripts/
# OCTOPUS_TOP=~/Octopus/octopus/

rm *.txt
octopus > log

$HELPER_DIR/extract_generic.sh static/info 'Eigenvalues'  'Dipole' | head -n -1 > info.txt

cp *.txt $OCTOPUS_TOP/doc/tutorials/octopus_basics/total_energy_convergence/2.methane/
