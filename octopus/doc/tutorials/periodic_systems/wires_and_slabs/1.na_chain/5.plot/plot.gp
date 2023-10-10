set xlabel "k_x/G_x"
set ylabel "E (eV)"
set t postscript enhanced color font "Monospace-Bold,25" landscape size 11,8.5
set output "Na_chain_bands_1D_3D.eps"
set key left

set rmargin 4.5
set lmargin 12.5
set tmargin 3.2
set bmargin 5.5

#set ytics 0.2
set mytics 2
set yrange [:6]
set zeroaxis

plot for [col=5:5+4] 'bandstructure_1D' u 2:(column(col)+3.06851682) w lp notitle lc 1 lt 7, for [col=5:5+4] 'bandstructure_3D' u 2:(column(col)+2.83078123) w lp notitle lc 2 lt 6, 20 w l t "1D-periodic" lc 1, 20 w l t "3D-periodic" lc 82
