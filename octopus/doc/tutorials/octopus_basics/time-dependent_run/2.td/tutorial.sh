#!/usr/bin/env bash

# these variables need to be defined:
# HELPER_DIR=~/HUGO/octopus-documentation/scripts/
# OCTOPUS_TOP=~/Octopus/octopus/

cp $OCTOPUS_TOP/testsuite/tutorials/07-octopus_basics-time_dependent_propagation.02-td.inp inp
octopus > log

$HELPER_DIR/extract_generic.sh log 'IonsConstantVelocity' 'Info:' | head -n -1 > input-params.txt
$HELPER_DIR/extract_generic.sh log ' Time-Dependent Simulation ' 'Info:' | head -n 7  > td-start.txt
$HELPER_DIR/extract_generic.sh log ' Time-Dependent Simulation ' 'Info:' | tail -n 5 | head -n -2  > td-end.txt

cp tutorial.sh *.txt log $OCTOPUS_TOP/doc/tutorials/octopus_basics/time-dependent_run/2.td/
