set xlabel "x (a.u.)"
set t postscript enhanced color font "Monospace-Bold,25" landscape size 11,8.5
set output "wavefunctions.eps"
set key top right

set rmargin 4.5
set lmargin 10.5
set tmargin 3.2
set bmargin 5.5

set yrange [-1:4]


p 'static/v0.y=0,z=0' w l ls 1 lw 3 t 'potential', \
  'static/wf-k001-st0001.y=0,z=0' w l ls 2 lw 3 t '1st state', \
  'static/wf-k001-st0002.y=0,z=0' w l ls 3 lw 3 t '2nd state'
