#!/bin/bash
list="10 15 20 25 30 35"
export OCT_PARSE_ENV=1

export OCT_ExtraStates=0
export OCT_FromScratch=yes
export OCT_CalculationMode=gs
$OCTPATH/octopus >& out-gs
unset OCT_FromScratch
for ExtraStates in $list
do
    export OCT_ExtraStates=$ExtraStates
    export OCT_CalculationMode=unocc
    $OCTPATH/octopus >& out-unocc-$ExtraStates
    export OCT_CalculationMode=casida
    $OCTPATH/octopus >& out-casida-$ExtraStates
    $OCTPATH/oct-casida_spectrum >& out-spec-$ExtraStates
    mv casida/spectrum.casida spectrum.casida-$ExtraStates
done
unset OCT_ExtraStates OCT_CalculationMode


