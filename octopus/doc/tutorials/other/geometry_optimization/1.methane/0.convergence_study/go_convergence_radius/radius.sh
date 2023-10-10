#!/bin/bash
echo "#Sp    Energy" > radius.log
list="4.0 5.0 6.0 7.0 8.0"
export OCT_PARSE_ENV=1
for Radius in $list
do
    export OCT_Radius=$(echo $Radius*1.8897261328856432 | bc)
    $OCTPATH/octopus >& out-$Radius
    energy=`grep Total static/info  | head -1 | cut -d "=" -f 2`
    echo $Radius $energy >> radius.log
    cp static/forces forces-$Radius
    cp min.xyz min.xyz-$Radius
    cp -r geom geom-$Radius
    rm -rf restart
done
 unset OCT_Radius
