#!/usr/bin/env bash

# these variables need to be defined:
# HELPER_DIR=~/HUGO/octopus-documentation/scripts/
# OCTOPUS_TOP=~/Octopus/octopus/

rm *.log *.txt

cp $OCTOPUS_TOP/testsuite/tutorials/07-octopus_basics-time_dependent_propagation.01-gs.inp inp


octopus > log


cp tutorial.sh  $OCTOPUS_TOP/doc/tutorials/octopus_basics/time-dependent_run/1.gs
