#!/bin/bash
echo "#CH distance   Total Energy" > ch.log
list="0.90 0.95 1.00 1.05 1.10 1.15 1.20 1.25 1.30"
for CH in $list
do
    sed -ie "/CH = /s|[0-9.]\+|$CH|" inp
    $OCTPATH/octopus >& out-$CH
    energy=`grep Total static/info  | head -1 | cut -d "=" -f 2`
    echo $CH $energy >> ch.log
    rm -rf restart
done

