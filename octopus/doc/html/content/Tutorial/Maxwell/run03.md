---
Title: "External currents and PML"
series: "Tutorials"
tutorials: "Maxwell"
author: ["Franco Bonafé","René Jestaedt","Heiko Appel"]
Weight: 14
---

### External current density

#### Not absorbing boundaries

Instead of an incoming external plane wave, Octopus can simulate also external
current densities placed inside the simulation box. In this example we place an
array of such current densities in the simulation box and study the properties
of the emitted field and the effect of the absorbing boundaries.

Since we start with no absorbing boundaries, we set the box size to 10.0.

{{< code-block >}}
#include_input doc/tutorials/maxwell/3.external-current/1.gaussian_current_pulse/inp
{{< /code-block >}}

The external current density is switched on by the corresponding options and
two blocks define its spatial distribution and its temporal behavior. In this
example we place three sources, located near the box boundary in the negative x
direction, and separated by 5 bohr along the y axis. As all of them are
polarized in the z direction, only the z component is non zero. The spatial
distribution of our example external current sources is a Gaussian distribution
in 3D.  The temporal pulse is a sinusoidal wave with a wavelength of 5 bohr,
modulated by a Gaussian envelope centered at half the total simulation time,
and with a variance $\sigma^2$ of 0.25.

{{< code-block >}}
#include_input_snippet doc/tutorials/maxwell/3.external-current/1.gaussian_current_pulse/inp current
{{< /code-block >}}

In addition to the electric field and energy density, we set the external
current as output, to visualize the shape and temporal profile and check that
our input file is correct. As we want a more frequent output for this quantity,
and only the z=0 plane output is sufficient, we select a specific output
format and interval for this quantity only. For the other quantities, the
global formats and interval apply.

{{< code-block >}}
#include_input_snippet doc/tutorials/maxwell/3.external-current/1.gaussian_current_pulse/inp output
{{< /code-block >}}

After running Octopus, we first visualize the current pulse using the following gnuplot script.

{{% expand "gnuplot script" %}}
```
set pm3d
set view map
set palette defined (-0.005 "blue", 0 "white", 0.005"red")
set term png size 1000,500

unset surface
unset key

set output 'plot1.png'

set origin 0.025,0
set size 0.3,0.9
set size square
set title 'External current density J_z - step 30'
sp [-10:10][-10:10] 'Maxwell/output_iter/td.0000030/external_current-z.z=0' u 1:2:3

set origin 0.35,0
set size 0.3,0.9
set size square
set title 'External current density J_z - step 90'
sp [-10:10][-10:10] 'Maxwell/output_iter/td.0000090/external_current-z.z=0' u 1:2:3

set origin 0.675,0
set size 0.3,0.9
set size square
set title 'External current density J_z - step 150'
sp [-10:10][-10:10] 'Maxwell/output_iter/td.0000150/external_current-z.z=0' u 1:2:3
unset multiplot

unset multiplot
```
{{% /expand %}}

{{< figure src="/images/Maxwell/tutorial_03.1-plot_current_2d.png" width="80%" >}}

We observe the three sources with their Gaussian spatial profile. Now, to see
the temporal profile, we extract the value of the current density at one
specific point, namely (-8 Bohr, 0, 0), from all of the output files.
Inspecting one of the output_iter/td.*/external_current-z.z\=0 files, we find
that the aforementioned point corresponds to line 195, and the field value
is in column 3.  Then, to grab these numbers for all time steps, we can run
this simple bash command:

{{< code-block >}}
for td in Maxwell/output_iter/td.0000*; do jj=$(sed -n '195p' $td//external_current-z.z\=0 | awk {'print $3'}); tdi=${td:27:6}; echo $tdi $jj >> current_vs_t.dat; done
{{< /code-block >}}

This creates a file with the time step number in the first column, and the
value of the current in the second column. In order to plot this, we can use
gnuplot:

{{% expand "gnuplot script" %}}
```
set output 'plot_j_td.png
set term png size 500,400
set xlabel 'time step'
set ylabel 'current density J_z'
p 'current_vs_t.dat' w l
```
{{% /expand %}}

{{< figure src="/images/Maxwell/tutorial_03.1-plot_current_td.png" width="50%" >}}

The EM irradiated by these current pulses should look like plane waves far
enough from the sources, and should have the same temporal profile than the
source. Let's look at the total electric field at the origin of the box, that
is printed to Maxwell/td.general/total_e_field_z:

{{% expand "gnuplot script" %}}
```
set output 'plot_e_field_at_orig.png
set term png size 500,400
set xlabel 'time step'
set ylabel 'Total E field E_z'
p 'Maxwell/td.general/total_e_field_z' u 1:3 w l
```
{{% /expand %}}

{{< figure src="/images/Maxwell/tutorial_03.1-plot_e_field_at_orig.png" width="50%" >}}

Clearly the amplitude increases and does not have a Gaussian envelope. Let's
inspect the E field in the z=0 plane to understand why, using the following
plotting script:

{{% expand "gnuplot script" %}}
```
set pm3d
set view map
set palette defined (-0.05 "blue", 0 "white", 0.05"red")
set term png size 1000,300

unset surface
unset key

set output 'plot_efield.png'

set xlabel 'x-direction'
set ylabel 'y-direction'
set cbrange [-0.01:0.01]

set multiplot

set origin 0.025,0
set size 0.3,0.9
set size square
set title 'Electric field E_z - step 60'
sp [-10:10][-10:10] 'Maxwell/output_iter/td.0000060/e_field-z.z=0' u 1:2:3

set origin 0.35,0
set size 0.3,0.9
set size square
set title 'Electric field E_z - step 90'
sp [-10:10][-10:10] 'Maxwell/output_iter/td.0000090/e_field-z.z=0' u 1:2:3

set origin 0.675,0
set size 0.3,0.9
set size square
set title 'Electric field E_z - step 120'
sp [-10:10][-10:10] 'Maxwell/output_iter/td.0000120/e_field-z.z=0' u 1:2:3
unset multiplot
```
{{% /expand %}}

{{< figure src="/images/Maxwell/tutorial_03.1-plot_e_field_2d.png" width="80%" >}}

After step 90 we start seeing a strong effect from the reflections that occur
at the box boundaries. In the following, we activate the PML boundary
conditions to improve this behaviour.

#### PML absorbing boundaries

Taking the previous input file as template, we now modify the absorbing
boundaries block:

{{< code-block >}}
#include_input_snippet doc/tutorials/maxwell/3.external-current/2.gaussian_current_pulse_with_pml/inp boundaries
{{< /code-block >}}

We set a boundary width of 5 Bohr (it should be at least 8 grid points, which
according to our spacing, would be 4 Bohr). With this additional absorbing
width, we have to update the simulation box dimensions.

{{< code-block >}}
#include_input_snippet doc/tutorials/maxwell/3.external-current/2.gaussian_current_pulse_with_pml/inp box
{{< /code-block >}}

The rest of the input options remain the same. After running the code, we now check the temporal profile
of the E field at the box center, plotting again column 3 vs. column 1 of the file 
Maxwell/td.general/total_e_field_z using the same script as before:

{{< figure src="/images/Maxwell/tutorial_03.2-plot_e_field_at_orig.png" width="50%" >}}

Now we do observe the Gaussian envelope. Looking at the z=0 plane E field, (plotting it
with the same script as before) we confirm the reflection problem is gone:

{{< figure src="/images/Maxwell/tutorial_03.2-plot_e_field_2d.png" width="80%" >}}

#### Assignment:

1. Is the E field only z-polarized? Check the other polarization directions.
2. In which direction is the B field predominantly polarized?
3. How much does the generated E_z field differ from a plane wave in the axis x = 8 Bohr? Grab the
necessary data from the output files.

{{< tutorial-footer >}}
