#!/bin/bash
echo "#Sp    Energy        s_eigen   p_eigen" > spacing.log
list="0.20 0.16 0.12 0.08 0.04"
export OCT_PARSE_ENV=1
for Spacing in $list
do
    export OCT_Spacing=$(echo $Spacing*1.8897261328856432 | bc)
    echo $Spacing
    octopus >& out-$Spacing
    energy=`grep Total static/info  | head -1 | cut -d "=" -f 2`
    echo $Spacing $energy $seigen $peigen >> spacing.log
    rm -rf restart
done
unset OCT_Spacing


