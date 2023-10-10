#!/usr/bin/env bash

# these variables need to be defined:
# HELPER_DIR=~/HUGO/octopus-documentation/scripts/
# OCTOPUS_TOP=~/Octopus/octopus/

octopus > log

./extract.pl > spectrum.txt

cp tutorial.sh extract.pl spectrum.txt  $OCTOPUS_TOP/doc/tutorials/optical_response/optical_spectra_from_sternheimer/2.em_resp
