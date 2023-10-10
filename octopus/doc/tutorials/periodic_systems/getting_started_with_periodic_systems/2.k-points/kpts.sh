#!/bin/bash
echo "#nk Total Energy" > kpts.log
list="2 4 6 8 10 12"
for nkpt in $list
do
    echo "$nkpt"
    sed -i "s/nk = [0-9]\+/nk = $nkpt/" inp
    octopus >& out-$nkpt
    energy=`grep Total static/info  | head -2 | tail -1 | cut -d "=" -f 2`
    echo $nkpt $energy >> kpts.log
    rm -rf restart
done
