set xlabel "Energy (eV)"
set ylabel "Strength Function (1/eV)"
set t postscript enhanced color font "Monospace-Bold,25" landscape size 11,8.5
set output "Absorption_spectrum_CH4_radius.eps"
set key top center title "Radius (angstrom)"

set rmargin 4.5
set lmargin 10.5
set tmargin 3.2
set bmargin 5.5

set yrange [0:3]
set xrange [7:15]
#set mxtics 2

plot "cross_section_vector-3.5" u 1:5 w l lw 1 t "3.5", "cross_section_vector-4.5" u 1:5 w l lw 1 t "4.5", "cross_section_vector-5.5" u 1:5 w l lw 1 t "5.5", "cross_section_vector-6.5" u 1:5 w l lw 1 t "6.5", "cross_section_vector-7.5" u 1:5 w l lw 1 t "7.5"
