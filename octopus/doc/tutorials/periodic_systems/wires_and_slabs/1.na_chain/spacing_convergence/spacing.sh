#!/bin/bash
echo "#Sp    Energy" > spacing.log
list="0.34 0.32 0.30 0.28 0.26 0.24"
export OCT_PARSE_ENV=1
for Spacing in $list
do
    export OCT_Spacing=$(echo $Spacing*1.8897261328856432 | bc)
    $OCTPATH/octopus >& out-$Spacing
    energy=`grep Total static/info  | head -2 | tail -1 | cut -d "=" -f 2`
    echo $Spacing $energy >> spacing.log
    rm -rf restart
done
 unset OCT_Spacing
