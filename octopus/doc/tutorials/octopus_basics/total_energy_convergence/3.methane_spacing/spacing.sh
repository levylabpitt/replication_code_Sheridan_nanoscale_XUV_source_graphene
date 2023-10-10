#!/bin/bash
echo "#Sp    Energy" > spacing.log
list="0.22 0.20 0.18 0.16 0.14 0.12 0.10"
export OCT_PARSE_ENV=1
for Spacing in $list
do
    export OCT_Spacing=$(echo $Spacing*1.8897261328856432 | bc)
    octopus >& out-$Spacing
    energy=`grep Total static/info  | head -1 | cut -d "=" -f 2`
    echo $Spacing $energy >> spacing.log
    rm -rf restart
done
 unset OCT_Spacing
