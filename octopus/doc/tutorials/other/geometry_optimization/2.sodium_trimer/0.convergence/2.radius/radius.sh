#!/bin/bash
echo "#Sp    Energy  Max. Force" > radius.log
list="5.0 6.0 7.0 8.0 9.0"
export OCT_PARSE_ENV=1
for Radius in $list
do
    export OCT_Radius=$(echo $Radius*1.8897261328856432 | bc)
    $OCTPATH/octopus >& out-$Radius
    energy=`grep Total static/info  | head -1 | cut -d "=" -f 2`
    max_force=`grep "Max" static/info | awk '{print $4, $5, $6}'`
    echo $Radius $energy $max_force >> radius.log
    rm -rf restart
done
 unset OCT_Radius
