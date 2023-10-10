#!/usr/bin/env bash

# these variables need to be defined:
# HELPER_DIR=~/HUGO/octopus-documentation/scripts/
# OCTOPUS_TOP=~/Octopus/octopus/

octopus > log

head -n 9 casida/casida > casida.txt 
head -n 9 casida/casida_excitations/00001 > casida-excitations.txt

cp tutorial.sh *.txt $OCTOPUS_TOP/doc/tutorials/optical_response/optical_spectra_from_casida/3.casida/
