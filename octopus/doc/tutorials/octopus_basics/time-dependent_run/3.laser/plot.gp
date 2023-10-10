set xlabel "t (hbar/eV)"
set ylabel "Amplitude (eV/Angstroms)"
set t postscript enhanced color font "Monospace-Bold,25" landscape size 11,8.5
set output "Tutorial_TD_Laser.eps"
unset key

set format y "%.1f"

set rmargin 4.5
set lmargin 10.5
set tmargin 3.2
set bmargin 5.5


plot "td.general/laser" u 2:3 w l lw 3
