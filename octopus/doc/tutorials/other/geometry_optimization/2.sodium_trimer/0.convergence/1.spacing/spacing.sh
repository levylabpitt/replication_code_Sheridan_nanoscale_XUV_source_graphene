#!/bin/bash
echo "#Sp    Energy  Max. Force" > spacing.log
list="0.32 0.30 0.28 0.26 0.24 0.22 0.20"
export OCT_PARSE_ENV=1
for Spacing in $list
do
    export OCT_Spacing=$(echo $Spacing*1.8897261328856432 | bc)
    $OCTPATH/octopus >& out-$Spacing
    energy=`grep Total static/info  | head -1 | cut -d "=" -f 2`
    max_force=`grep "Max" static/info | awk '{print $4, $5, $6}'`
    echo $Spacing $energy $max_force >> spacing.log
    rm -rf restart
done
 unset OCT_Spacing
