set xlabel "Spacing (Angstroms)"
set ylabel "Energy (eV)"
set t postscript enhanced color font "Monospace-Bold,25" landscape size 11,8.5
set output "NitrogenSpacing.eps"
set key left

set rmargin 4.5
set lmargin 9.5
set tmargin 3.2
set bmargin 5.5


plot "spacing.log" u 1:($2+261.78536939) w lp ps 2 pt 6 lw 3 t "total energy", "" u 1:($3+18.389733) w lp ps 2 pt 7 lw 2 t "s-eigen", "" u 1:($4+7.248998) w lp ps 2 pt 6 lw 2 t "p-eigen"
