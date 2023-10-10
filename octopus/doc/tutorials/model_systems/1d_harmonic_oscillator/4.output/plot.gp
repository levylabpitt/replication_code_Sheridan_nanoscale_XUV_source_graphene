set xlabel "x (a.u.)"
set t postscript enhanced color font "Monospace-Bold,25" landscape size 11,8.5
set output "HO_wavefunctions.eps"
set key top left

set rmargin 4.5
set lmargin 10.5
set tmargin 3.2
set bmargin 5.5


plot 'static/wf-st0001.y=0,z=0' u 1:(abs($2)) t 'n=1' w l, \
     'static/wf-st0002.y=0,z=0' t 'n=2' w l, \
     'static/wf-st0004.y=0,z=0' t 'n=4' w l, \
     'static/wf-st0006.y=0,z=0' t 'n=6' w l
