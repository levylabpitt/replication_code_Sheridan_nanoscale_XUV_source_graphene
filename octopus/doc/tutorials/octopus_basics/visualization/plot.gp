set xlabel "x"
set ylabel "y"
set zlabel "density"

set t postscript enhanced color font "Monospace-Bold,20" landscape size 11,8.5
set output "benzene_2D.eps"
set key left

set rmargin 4.5
set lmargin 9.5
set tmargin 3.2
set bmargin 5.5

set view 20, 350

splot [-5:5][-5:5][0:5] "static/density.z=0" w pm3d t ''
