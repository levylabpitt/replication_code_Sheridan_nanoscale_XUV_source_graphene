unset xlabel
set ylabel "E (hartree)"
set t postscript enhanced color font "Monospace-Bold,25" landscape size 11,8.5
set output "Si_bandstructure.eps"
unset key

set rmargin 4.5
set lmargin 12.5
set tmargin 3.2
set bmargin 5.5

set ytics 0.2
set mytics 2

set xtics ("{/Helvetica L}" 0.0, "{/Symbol G}" 0.26401290, "{/Helvetica X}" 0.56886874, "{/Symbol G}" 1.0)

plot for [col=5:5+4] 'static/bandstructure' u 1:(column(col)-0.149644) w l notitle ls 1, for [col=5+4:5+9] 'static/bandstructure' u 1:(column(col)-0.149644) w l notitle ls 2
