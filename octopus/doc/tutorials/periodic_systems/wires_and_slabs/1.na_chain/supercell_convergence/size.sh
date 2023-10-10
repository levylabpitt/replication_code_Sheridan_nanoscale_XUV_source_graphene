#!/bin/bash
echo "#Size Total Energy" > size.log
list="5.29 6.0 7.0 8.0 9.0"
for size in $list
do
    sed -i "s/Size = [0-9.]\+/Size = $size/" inp
    $OCTPATH/octopus >& out-$size
    energy=`grep Total static/info  | head -2 | tail -1 | cut -d "=" -f 2`
    echo $size $energy >> size.log
    rm -rf restart
done
