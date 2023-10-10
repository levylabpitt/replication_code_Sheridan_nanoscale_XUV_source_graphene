#!/usr/bin/env bash

# these variables need to be defined:
# HELPER_DIR=~/HUGO/octopus-documentation/scripts/
# OCTOPUS_TOP=~/Octopus/octopus/

cp inp_x inp
oct-propagation_spectrum > spectrum.log 2> spectrum.err

$HELPER_DIR/extract_generic.sh spectrum.log ' Spectrum Options ' 'Octopus emitted'  > spectrum-out.txt
$HELPER_DIR/extract_generic.sh spectrum.err 'Parser warning' '                '  > parser-warnings.txt
head -n 21 td.general/multipoles.1 > multipoles-head.txt

$HELPER_DIR/extract_generic.sh cross_section_vector.1 '# nspin ' '#%'  > cross_section_vector-1.txt
$HELPER_DIR/extract_generic.sh cross_section_vector.1 '# Number ' '#%'  > cross_section_vector-2.txt
$HELPER_DIR/extract_generic.sh cross_section_vector.1 '# Electronic ' '#%'  > cross_section_vector-3.txt
$HELPER_DIR/extract_generic.sh cross_section_vector.1 '    Energy    ' '#%' | head -n 5 > cross_section_vector-4.txt

head -n 17 cross_section_tensor > cross_section_tensor-head.txt

cp tutorial.sh *.txt log $OCTOPUS_TOP/doc/tutorials/optical_response/optical_spectra_from_time-propagation/3.propagation_spectra/
