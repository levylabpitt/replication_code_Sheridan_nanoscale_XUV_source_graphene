#!/usr/bin/env bash

# these variables need to be defined:
# HELPER_DIR=~/HUGO/octopus-documentation/scripts/
# OCTOPUS_TOP=~/Octopus/octopus/

octopus > log

cp static/eigenvalues eigenvalues.txt

cp tutorial.sh *.txt $OCTOPUS_TOP/doc/tutorials/optical_response/optical_spectra_from_casida/2.unocc/
