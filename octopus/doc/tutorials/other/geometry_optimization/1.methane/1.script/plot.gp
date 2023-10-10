set xlabel "CH-bond distance (Angstroms)"
set ylabel "Total Energy (eV)"
set t postscript enhanced color font "Monospace-Bold,25" landscape size 11,8.5
set output "Methane_Energy_vs_CH.eps"
unset key

set rmargin 4.5
set lmargin 12.5
set tmargin 3.2
set bmargin 5.5


plot "ch.log" u 1:2 w lp ps 2 pt 6 lw 3
