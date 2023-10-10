set ylabel "E (hartree)"
set t postscript enhanced color font "Monospace-Bold,25" landscape size 11,8.5
set output "Si_unfolded.eps"
unset key

set palette defined (0 "white", 100 "dark-violet", 600 "yellow")

set xtics ("{/Helvetica L}" 0.0, "{/Symbol G}" 0.26401290, "{/Helvetica X}" 0.56886874, "{/Symbol G}" 1.0)

plot 'static/ake.dat' u 1:2:3 w image  notitle
