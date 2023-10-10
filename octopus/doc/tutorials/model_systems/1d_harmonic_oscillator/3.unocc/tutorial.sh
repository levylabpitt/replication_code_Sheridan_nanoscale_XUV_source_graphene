#!/usr/bin/env bash

# these variables need to be defined:
# HELPER_DIR=~/HUGO/octopus-documentation/scripts/
# OCTOPUS_TOP=~/Octopus/octopus/

octopus > log

cp static/eigenvalues eigenvalues.txt



cp tutorial.sh eigenvalues.txt $OCTOPUS_TOP/doc/tutorials/model_systems/1d_harmonic_oscillator/3.unocc/
