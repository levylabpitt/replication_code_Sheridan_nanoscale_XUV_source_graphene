---
Title: "Plane waves in vacuum"
series: "Tutorials"
tutorials: "Maxwell"
author: ["Franco Bonafé","René Jestaedt","Heiko Appel"]
Weight: 11
---

## Cosinoidal plane wave in vacuum

In this tutorial, we will describe the propagation of waves in vacuum. Let's
start with one wave with a cosinoidal envelope.

{{< code-block >}}
#include_input doc/tutorials/maxwell/1.free-propagation/1.1_pulse_td/inp
{{< /code-block >}}

The total size of the (physical) simulation box is 20.0 x 20.0 x 20.0 Bohr (10
Bohr half length in each direction) with a spacing of 0.5 Bohr in each
direction. As discussed in the input overview, also the boundary points have to
be accounted for when defining the box size in the input file. In other words,
the user has to know how much space these points will add, which is given by
the derivatives order (in this case, 4) times the spacing. This is 4 * 0.5 bohr
= 2 bohr.  Hence the total size of the simulation box is chosen to be 12.0 bohr
in each direction.

{{< code-block >}}
#include_input_snippet doc/tutorials/maxwell/1.free-propagation/1.1_pulse_td/inp box
{{< /code-block >}}

The boundary conditions are chosen as plane waves to simulate the incoming wave
without any absorption at the boundaries.

{{< code-block >}}
#include_input_snippet doc/tutorials/maxwell/1.free-propagation/1.1_pulse_td/inp boundaries
{{< /code-block >}}

The simulation time step is set to the value given by the Courant condition.
For equal spacing in the three dimensions, this is equal to {{< code-inline >}} dt = spacing /
(sqrt(3)*c) = 0.002106 {{< /code-inline >}} atomic units of time, in this case,
c being the speed of light in atomic units (137.03599).
The system propagates 150 time steps = 0.316 atomic units.

{{< code-block >}}
#include_input_snippet doc/tutorials/maxwell/1.free-propagation/1.1_pulse_td/inp timestep
{{< /code-block >}}

The incident cosinoidal plane is evaluated at the boundaries to be fed into the
simulation box. The pulse has a width of 10.0 bohr, a spatial shift of 25.0
Bohr in the negative x-direction, an amplitude of 0.05 a.u., and a wavelength
of 10.0 bohr.

{{< code-block >}}
#include_input_snippet doc/tutorials/maxwell/1.free-propagation/1.1_pulse_td/inp field
{{< /code-block >}}

Finally, the output options are set:

{{< code-block >}}
#include_input_snippet doc/tutorials/maxwell/1.free-propagation/1.1_pulse_td/inp output
{{< /code-block >}}

Contour plot of the electric field in z-direction after 50 time steps (t=0.105)
and 100 time steps (t=0.21):

{{% expand "Example gnuplot script" %}}
```bash
set pm3d
set view map
set palette defined (-0.1 "blue", 0 "white", 0.1 "red")
set term png size 1000,500

unset surface
unset key

set output 'plot.png'

set xlabel 'x-direction'
set ylabel 'y-direction'
set cbrange [-0.1:0.1]

set multiplot

set origin 0.025,0
set size 0.45,0.9
set size square
set title 'Electric field E_z (t=0.105328 au)'
sp [-10:10][-10:10][-0.1:0.1] 'Maxwell/output_iter/td.0000050/e_field-z.z=0' u 1:2:3

set origin 0.525,0
set size 0.45,0.9
set size square
set title 'Electric field E_z (t=0.210656 au)'
sp [-10:10][-10:10][-0.1:0.1] 'Maxwell/output_iter/td.0000100/e_field-z.z=0' u 1:2:3

unset multiplot
```
If you copy this script into a file called {{< file "plot.gnu" >}}, you can create the plots by:
```bash
gnuplot plot.gnu
```
{{% /expand %}}

{{< figure src="/images/Maxwell/tutorial_01_run_electric_field_contour.png" width="50%" >}}


Maxwell fields at the origin and Maxwell energy inside the free Maxwell
propagation region of the simulation box:

{{% expand "gnuplot script" %}}
```bash
set term png size 1000,400
set output 'plot2.png'

set multiplot
set style data l

set xlabel 'time [a.u.]'

set origin 0.025,0
set size 0.45,0.9
set title 'Electric and magnetic fields'
set ylabel 'Maxwell fields [a.u.]'
p 'Maxwell/td.general/total_e_field_z' u 2:3 t 'E_z', 'Maxwell/td.general/total_b_field_y' u 2:($3*100) t '100 * B_y'

set origin 0.525,0
set size 0.45,0.9
set title 'Maxwell energy'
set ylabel 'Maxwell energy [a.u.]'
p  'Maxwell/td.general/maxwell_energy' u 2:3  t 'Maxwell energy'

unset multiplot

```
{{% /expand %}}

{{< figure src="/images/Maxwell/tutorial_01_run_maxwell_energy_and_fields.png" width="50%" >}}

## Interference of two cosinoidal plane waves

Instead of only one plane wave, we simulate in this example two different plane
waves with different wave-vectors entering the simulation box, interfering and
leaving the box again. In addition to the wave from the last tutorial, we add a
second wave with different wave length, and entering the box at an angle, and
shifted by 28 bohr along the corresponding direction of propagation.  Both
electric fields are polarized only in z-direction, and the magnetic field only
in y-direction.

{{< code-block >}}
#include_input doc/tutorials/maxwell/1.free-propagation/2.2_pulses_td/inp
{{< /code-block >}}

Using a similar scipt than in the previous example, the contour plots for the electric field in the z=0 plane can be obtained.

Contour plot of the electric field in z-direction after 50 time steps for
t=0.11 and 100 time steps for t=0.21:
{{< figure src="/images/Maxwell/tutorial_02_run_electric_field_contour.png" width="50%" >}}

Fields at the origin and total energy inside the propagation region of the simulation box:
{{< figure src="/images/Maxwell/tutorial_02_run_maxwell_energy_and_fields.png" width="50%" >}}


{{< tutorial-footer >}}
