---
Title: "Dispersive linear media"
series: "Tutorials"
tutorials: "Maxwell"
author: ["Franco Bonafé","Heiko Appel"]
Weight: 20
---

In addition to the possibility to include objects described as static linear
media, Octopus allows for simulations with dispersive media. For example, let's
consider the following input file, for propagation of a pulse through a sphere
described by a Drude polarizability:

{{< code-block >}}
#include_input doc/tutorials/maxwell/4.drude-medium/inp
{{< /code-block >}}

Here, we need to define the variables {{<code "MediumPoleEnergy">}} and
{{<code "MediumPoleDamping">}}, which represent the plasma frequency $\omega_p$ and
inverse lifetime $\gamma$ of a Drude pole. For this example, the parameters are those
for the Drude peak of gold, as taken from the literature [^footnote-1]. Also,
we can define {{<code "MediumCurrentCoordinates">}} to obtain the
polarization current at certain points. The rest of the input variables of
the medium are the same that have been used for the static linear medium. In
this case, we are using an off file that contains the shape of a sphere of 80
nm radius, which is displaced from the origin by one radius to the negative x
direction.

{{% expand "off file" %}}
{{< code-block >}}
#include_input doc/tutorials/maxwell/4.drude-medium/gold-np-r80nm.off
{{< /code-block >}}
{{% /expand %}}

For practical reasons, we set the box size as a function of the incoming pulse
wavelength {{<code "l_zero">}}, and we make the box larger in the direction of propagation.
We also add the spacing needed for the PML boundaries, defined below. As you
can see, this is defined as 550 nm, using the default parameter nm to convert
it to atomic units. Also, we set a relatively small Courant number, which will
replace the $1/\sqrt(3)$ value that was explained before. This is because for
Drude media, the time step is also limited by the plasma frequency of the
metal, and not only by the grid spacing, so a larger value will cause the
simulation to explode (you are welcome to try and increase it, plot the total
integrated energy, and observe the divergence after the threshold). For this
setup, 0.1 is an appropriate Courant number.

After we run the simulation, we can plot again the z-component of the electric
field in the xz-plane for three different time steps, using the following script:

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
set cbrange [-0.8:0.8]

set multiplot

set origin 0.025,0
set size 0.3,0.9
set size square
set title 'Electric field E_z - step 60'
sp [-400:400][-400:400] 'Maxwell/output_iter/td.0000060/e_field-z.y=0' u ($1/18.897):($2/18.897):3

set origin 0.35,0
set size 0.3,0.9
set size square
set title 'Electric field E_z - step 120'
sp [-400:400][-400:400] 'Maxwell/output_iter/td.0000120/e_field-z.y=0' u ($1/18.897):($2/18.897):3

set origin 0.675,0
set size 0.3,0.9
set size square
set title 'Electric field E_z - step 240'
sp [-400:400][-400:400] 'Maxwell/output_iter/td.0000240/e_field-z.y=0' u ($1/18.897):($2/18.897):3
unset multiplot
```
{{% /expand %}}

The plot we get is the following:

{{< figure src="/images/Maxwell/tutorial_04.1-drude-01.png" width="80%" >}}

As it is expected, as the medium reacts to the external field via its
polarization current density, screening it inside the sphere, as it is expected
from a metal. We can also note that the space discretization used is a rather
coarse mesh in this tutorial, due to computation time, but the outcomes are
still reasonable.

Finally, we can examine how these currents arise in time. Using the following
script we can plot the current at three different points (the ones we requested
in the input file, namely: near the surface of the sphere towards the negative
x axis, in the middle, and near the surface on the positive x axis). Also we
plot the E field at these points.

{{% expand "gnuplot script" %}}
```
set term png size 600,600

set output 'plot_j.png'
set xlabel 'time step'
set ylabel 'J_z'
set size square
plot 'NP/td.general/current_at_points.dat' u 1:5 w l title "before", 'NP/td.general/current_at_points.dat' u 1:8 w l title "middle", 'NP/td.general/current_at_points.dat' u 1:11 w l title "after"

set output 'plot_e.png'
set xlabel 'time step'
set ylabel 'E_z'
set size square
plot 'Maxwell/td.general/total_e_field_z' u 1:3 w l title "before", 'Maxwell/td.general/total_e_field_z' u 1:4 w l title "middle", 'Maxwell/td.general/total_e_field_z' u 1:5 w l title "after"
```
{{% /expand %}}

{{< figure src="/images/Maxwell/tutorial_04.1-drude-02.png" width="30%" >}}
{{< figure src="/images/Maxwell/tutorial_04.1-drude-03.png" width="30%" >}}

As can be seen, the current in the direction of the incident electric field
(called "before") is larger, screening the field effectively. This can be seen
by plotting the electric field in the other points, where the values in the
"middle" and "after" are much smaller than "before" (actually, the plotted
field has already been quenched by the medium, otherwise the waveform would be
that of the Gaussian envelope used). In addition to scattering, there is
absorption of energy by the nanoparticle, given by the imaginary part of the
polarizability. Using a setup like this, and processing the proper values of
the EM field in full space and time, it would be possible to calculate an
extinction spectrum of the sphere, which can be compared to the one calculated
using Mie theory, or other methods.


---------------------------------------------
[^footnote-1]: {{< article title="Optical properties of metallic films for vertical-cavity optoelectronic devices" authors="Aleksandar D. Rakić, Aleksandra B. Djurišić, Jovan M. Elazar, and Marian L. Majewski" journal="Applied Optics" volume="37" pages="5271" year="1998" url="https://opg.optica.org/ao/fulltext.cfm?uri=ao-37-22-5271&id=61190" doi="10.1364/AO.37.005271" link="AO" >}}

{{< tutorial-footer >}}
