set xlabel "k [(a+b)/pi]"
set ylabel "Energy [Hartree]"
set t postscript enhanced color font "Monospace-Bold,25" landscape size 11,8.5
set output "bandstructure.eps"
set key top left

set rmargin 4.5
set lmargin 10.5
set tmargin 3.2
set bmargin 5.5

plot for [col=5:5+9] 'static/bandstructure' u 1:(column(col)) w l notitle ls 1
