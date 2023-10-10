#!/usr/bin/env bash

# these variables need to be defined:
# HELPER_DIR=~/HUGO/octopus-documentation/scripts/
# OCTOPUS_TOP=~/Octopus/octopus/

oct-casida_spectrum > casida.log

cp casida/spectrum.casida spectrum.casida-singlets

head -n 9 casida/casida > casida.txt 
head -n 9 casida/casida_excitations/00001 > casida-excitations.txt

gnuplot plot.gp

cp spectrum.casida-singlets plot.gp tutorial.sh *.txt Absorption_spectrum_CH4_casida.eps $OCTOPUS_TOP/doc/tutorials/optical_response/optical_spectra_from_casida/4.spectrum
