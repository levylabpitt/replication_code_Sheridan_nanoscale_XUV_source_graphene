set xlabel "Energy (eV)"
set ylabel "Absorption"
set t postscript enhanced color font "Monospace-Bold,25" landscape size 11,8.5
set output "Eps2_Si.eps"
unset key

set rmargin 4.5
set lmargin 10.5
set tmargin 3.2
set bmargin 5.5

set xzeroaxis
set xrange [0:15]
set yrange [-20:]
set mxtics 5
set mytics 5
plot 'td.general/dielectric_function' u ($1*27.2114):3 w l lw 3

