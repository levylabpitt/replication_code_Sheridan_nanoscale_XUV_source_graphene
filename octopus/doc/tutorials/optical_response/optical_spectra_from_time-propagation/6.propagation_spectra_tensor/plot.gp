set xlabel "Energy (eV)"
set ylabel "Strength Function (1/eV)"
set t postscript enhanced color font "Monospace-Bold,25" landscape size 11,8.5
set output "Absorption_spectrum_CH4.eps"
unset key

set rmargin 4.5
set lmargin 10.5
set tmargin 3.2
set bmargin 5.5

set yrange [0:3]
set mxtics 5
plot "cross_section_vector.1" u 1:5 w l lw 3
