#!/bin/bash
list="0.26 0.24 0.22 0.20 0.18"
export OCT_PARSE_ENV=1
for Spacing in $list
do
    export OCT_Spacing=$(echo $Spacing*1.8897261328856432 | bc)
    export OCT_CalculationMode=gs
    octopus >& out-gs-$Spacing
    export OCT_CalculationMode=td
    octopus >& out-td-$Spacing
    oct-propagation_spectrum >& out-spec-$Spacing
    mv cross_section_vector cross_section_vector-$Spacing
    rm -rf restart
done
unset OCT_Spacing OCT_CalculationMode


