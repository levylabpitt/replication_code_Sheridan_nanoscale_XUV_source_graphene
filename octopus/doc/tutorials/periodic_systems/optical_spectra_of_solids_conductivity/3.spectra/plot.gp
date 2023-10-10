set xlabel "Energy (eV)"
set ylabel "Absorption"
set t postscript enhanced color font "Monospace-Bold,25" landscape size 11,8.5
set output "Si_abs_conduc.eps"
unset key

set rmargin 4.5
set lmargin 10.5
set tmargin 3.2
set bmargin 5.5

set xzeroaxis
set xrange [0:15]
set yrange [-20:200]
set mxtics 5
set mytics 5
plot 'td.general/conductivity' u ($1*27.2114):($2/$1*4*3.14159265*(10.18**3/4)/0.01) w l lw 3

