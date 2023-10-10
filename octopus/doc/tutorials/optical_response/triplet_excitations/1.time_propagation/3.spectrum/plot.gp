set xlabel "Energy (eV)"
set ylabel "Strength Function (1/eV)"
set t postscript enhanced color font "Monospace-Bold,25" landscape size 11,8.5
set output "Singlet_triplet_spectrum_CH4.eps"
set key top left

set rmargin 4.5
set lmargin 10.5
set tmargin 3.2
set bmargin 5.5

set yrange [0:]
set xrange [0:20]
set ytics 0.5
set mxtics 5
set mytics 5

plot "cross_section_vector-singlets" u 1:5 w l lw 2 t "singlet", "cross_section_vector-triplets" u 1:($8-$9) w l lw 2 t "triplet"
