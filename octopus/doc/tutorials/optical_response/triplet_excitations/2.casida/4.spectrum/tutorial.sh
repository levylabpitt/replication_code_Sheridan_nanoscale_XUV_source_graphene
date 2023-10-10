#!/usr/bin/env bash

# these variables need to be defined:
# HELPER_DIR=~/HUGO/octopus-documentation/scripts/
# OCTOPUS_TOP=~/Octopus/octopus/

oct-casida_spectrum > spectrum.log 

cp casida/spectrum.casida spectrum.casida-triplets
cp $OCTOPUS_TOP/doc/tutorials/optical_response/optical_spectra_from_casida/4.spectrum/spectrum.casida-singlets .

gnuplot plot.gp

cp tutorial.sh plot.gp Singlet_triplet_spectrum_Casida_CH4.eps $OCTOPUS_TOP/doc/tutorials/optical_response/triplet_excitations/2.casida/4.spectrum/
