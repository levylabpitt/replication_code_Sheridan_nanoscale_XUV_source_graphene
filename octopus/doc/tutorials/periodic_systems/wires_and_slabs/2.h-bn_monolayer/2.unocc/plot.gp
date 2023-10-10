unset xlabel
set ylabel "E (eV)"
set t postscript enhanced color font "Monospace-Bold,25" landscape size 11,8.5
set output "Tutorial_band_structure_HBN.eps"
unset key

set rmargin 4.5
set lmargin 12.5
set tmargin 3.2
set bmargin 5.5

set mytics 2

set xtics ("{/Symbol G}" 0.0, "{/Helvetica K}" 0.42264973, "{/Helvetica M}" 0.63397460, "{/Symbol G}" 1.0)

plot for [col=5:5+3] 'static/bandstructure' u 1:(column(col)+5.97340263) w l notitle ls 1, for [col=5+4:5+8] 'static/bandstructure' u 1:(column(col)+5.97340263) w l notitle ls 2
