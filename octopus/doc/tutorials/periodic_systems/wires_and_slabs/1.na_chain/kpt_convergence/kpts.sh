#!/bin/bash
echo "#nk Total Energy" > kpts.log
list="3 5 7 9 11 13 15"
for nkpt in $list
do
    sed -i "s/nk = [0-9]\+/nk = $nkpt/" inp
    $OCTPATH/octopus >& out-$nkpt
    energy=`grep Total static/info  | head -2 | tail -1 | cut -d "=" -f 2`
    echo $nkpt $energy >> kpts.log
    rm -rf restart
done
