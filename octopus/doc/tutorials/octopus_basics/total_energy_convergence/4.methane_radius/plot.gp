set xlabel "Radius (Angstroms)"
set ylabel "Total Energy (eV)"
set t postscript enhanced color font "Monospace-Bold,25" landscape size 11,8.5
set output "Radius_CH4.eps"
unset key

set format y "%.2f"

set rmargin 4.5
set lmargin 12.5
set tmargin 3.2
set bmargin 5.5


plot "radius.log" u 1:2 w lp ps 2 pt 6 lw 3
