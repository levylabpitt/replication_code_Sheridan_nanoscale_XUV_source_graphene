set xlabel "Energy (eV)"
set ylabel "Strength Function (1/eV)"
set t postscript enhanced color font "Monospace-Bold,25" landscape size 11,8.5
set output "Absorption_spectrum_CH4_spacing.eps"
set key top center title "Spacing (angstrom)"

set rmargin 4.5
set lmargin 10.5
set tmargin 3.2
set bmargin 5.5

set yrange [0:3]
set xrange [8.5:14.5]
set mxtics 2
#plot "cross_section_vector-0.18" u 1:5 w l lw 1 t "0.18", "cross_section_vector-0.20" u 1:5 w l lw 1 t "0.20", "cross_section_vector-0.22" u 1:5 w l lw 1 t "0.22", "cross_section_vector-0.24" u 1:5 w l lw 1 t "0.24", "cross_section_vector-0.26" u 1:5 w l lw 1 t "0.26"

plot "cross_section_vector-0.18" u 1:5 w l lw 1 t "0.18", "cross_section_vector-0.22" u 1:5 w l lw 1 t "0.22", "cross_section_vector-0.24" u 1:5 w l lw 1 t "0.24"
