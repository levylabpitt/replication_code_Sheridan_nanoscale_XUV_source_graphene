#!/bin/bash
list="8.5 9.5"
export OCT_PARSE_ENV=1
for Radius in $list
do
    export OCT_Radius=$(echo $Radius*1.8897261328856432 | bc)
    export OCT_CalculationMode=gs
    octopus >& out-gs-$Radius
    export OCT_CalculationMode=td
    octopus >& out-td-$Radius
    oct-propagation_spectrum >& out-spec-$Radius
    mv cross_section_vector cross_section_vector-$Radius
    rm -rf restart
done
unset OCT_Radius OCT_CalculationMode


