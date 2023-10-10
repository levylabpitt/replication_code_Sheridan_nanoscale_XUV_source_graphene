set xlabel "x1"
set ylabel "x2"
set t postscript enhanced color font "Monospace-Bold,20" landscape size 11,8.5
set key top left

set hidden3d
set pm3d
set contour
set ticslevel 0
unset key
unset surface

set output "Wfs-st0001.eps"
splot 'static/wf-st0001.z=0' using 1:2:3

set output "Wfs-st0002.eps"
splot 'static/wf-st0002.z=0' using 1:2:3

set output "Wfs-st0003.eps"
splot 'static/wf-st0003.z=0' using 1:2:3
