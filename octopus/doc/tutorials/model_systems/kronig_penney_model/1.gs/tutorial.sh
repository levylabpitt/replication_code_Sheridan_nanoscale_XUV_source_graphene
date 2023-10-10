#!/usr/bin/env bash

# these variables need to be defined:
# HELPER_DIR=~/HUGO/octopus-documentation/scripts/
# OCTOPUS_TOP=~/Octopus/octopus/

octopus > log

gnuplot plot.gp

cp inp tutorial.sh plot.gp *.eps  $OCTOPUS_TOP/doc/tutorials/model_systems/kronig_penney_model/1.gs/
