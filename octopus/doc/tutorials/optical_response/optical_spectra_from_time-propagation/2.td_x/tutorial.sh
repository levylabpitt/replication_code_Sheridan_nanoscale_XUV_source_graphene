#!/usr/bin/env bash

# these variables need to be defined:
# HELPER_DIR=~/HUGO/octopus-documentation/scripts/
# OCTOPUS_TOP=~/Octopus/octopus/

cp inp_x inp
echo "run x"
octopus > log

$HELPER_DIR/extract_generic.sh log 'Applying delta' 'kick mode:'  > kick-info.txt
head -n 21 td.general/multipoles > multipoles-head.txt

cp inp_? tutorial.sh *.txt log $OCTOPUS_TOP/doc/tutorials/optical_response/optical_spectra_from_time-propagation/2.td_x/

mv td.general/multipoles td.general/multipoles.1

cp inp_y inp
echo "run y"
octopus > log
mv td.general/multipoles td.general/multipoles.2


cp inp_z inp
echo "run z"
octopus > log
mv td.general/multipoles td.general/multipoles.3
