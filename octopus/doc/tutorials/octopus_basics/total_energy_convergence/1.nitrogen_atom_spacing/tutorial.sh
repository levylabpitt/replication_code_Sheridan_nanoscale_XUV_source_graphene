#!/usr/bin/env bash

# these variables need to be defined:
# HELPER_DIR=~/HUGO/octopus-documentation/scripts/
# OCTOPUS_TOP=~/Octopus/octopus/

cp $OCTOPUS_TOP/testsuite/tutorials/01-octopus_basics-getting_started.01-H_atom.inp inp

. ./spacing.sh

gnuplot plot.gp

cp *.eps tutorial.sh *.log $OCTOPUS_TOP/doc/tutorials/octopus_basics/total_energy_convergence/1.nitrogen_atom_spacing/
