#!/usr/bin/env bash

oct-propagation_spectrum > spectrum.log
head -n 20 cross_section_tensor > cross_section_tensor-head.txt

cp tutorial.sh *.txt $OCTOPUS_TOP/doc/tutorials/optical_response/use_of_symmetries_in_optical_spectra_from_time-propagation/3.spectrum/
