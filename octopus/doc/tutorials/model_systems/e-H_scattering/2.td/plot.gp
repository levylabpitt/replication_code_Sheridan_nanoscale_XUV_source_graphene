set terminal pngcairo font "Monospace-Bold,20" size 600,700

do for [t=0:1112:50] {

 set output sprintf('potential_%04i.png',t)

 set multiplot

 set xrange [-50:50]

 set tmargin 1.5
 set bmargin 1.5

 set size 1.0, 0.6
 set origin 0.0, 0.4
 set ylabel "n(x)" offset 2.0,0

 au2fs = 0.024189
 dt = 0.02973;

 set label sprintf('t=%.2f fs', t*dt*au2fs) at -15,0.45
 set yrange [0:0.5]
 plot sprintf('output_iter/td.%07i/density-sp1.y=0,z=0', t) u 1:2 w l lw 2 notitle
  
 unset label
 set size 1.0, 0.4
 set origin 0.0, 0.0
 set yrange [-0.5:0.0]
 set ylabel "v_{xc}(x)" offset 2.0,0
 set xlabel "x [Bohr]" offset 0, 0.6
 plot sprintf('output_iter/td.%07i/vxc-sp1.y=0,z=0', t) u 1:2 w l ls 2 lw 2 notitle

 unset label
 unset multiplot
}
