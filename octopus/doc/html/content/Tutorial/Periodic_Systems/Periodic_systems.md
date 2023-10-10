---
title: "Getting started with periodic systems"
series: "Tutorials"
tutorials: ["Periodic Systems"]
difficulties: "basic"
difficulties_weight: 1
theories: "DFT"
calculation_modes: ["Ground state", "Unoccupied"]
system_types: "bulk"
species_types: "Pseudopotentials"
features: ["Band structure", "DOS"]
description: "How to perform some basic calculation using bulk silicon as an example."
weight: 1
---


The extension of a ground-state calculation to a periodic system is quite straightforward in {{< octopus >}}. In this tutorial we will explain how to perform some basic calculation using bulk silicon as an example.

## Input
As always, we will start with a simple input file. In this case we will use a primitive cell of Si, composed of two atoms. 

{{< code-block >}}
#include_input testsuite/tutorials/06-octopus_basics-periodic_systems.01-silicon.inp
{{< /code-block >}}


Lets see more in detail some of the input variables:

* {{< code-inline >}}{{< variable "PeriodicDimensions" >}} = 3{{< /code-inline >}}: this input variable must be set equal to the number of dimensions you want to consider as periodic. Since the system is 3D ({{< code-inline >}}{{< variable "Dimensions" >}} = 3{{< /code-inline >}} is the default), by setting this variable to 3 we impose periodic boundary conditions at all borders. This means that we have a fully periodic infinite crystal.

* {{< code-inline >}}{{< variable "LatticeVectors" >}}{{< /code-inline >}} and {{< code-inline >}}{{< variable "LatticeParameters" >}}{{< /code-inline >}}: these two blocks are used to define the primitive lattice vectors that determine the unit cell. {{< code-inline >}}{{< variable "LatticeVectors" >}}{{< /code-inline >}} defines the direction of the vectors, while {{< variable "LatticeParameters" >}} defines their length.

* {{< code-inline >}}{{< variable "ReducedCoordinates" >}}{{< /code-inline >}}: the position of the atoms inside the unit cell, in reduced coordinates.

* {{< code-inline >}}{{< variable "KPointsGrid" >}}{{< /code-inline >}}: this specifies the ''k''-point grid to be used in the calculation. Here we employ a 2x2x2 Monkhorst-Pack grid with four shifts. The first line of the block defines the number of ''k''-points along each axis in the Brillouin zone. Since we want the same number of points along each direction, we have defined the auxiliary variable {{< code "nk = 2" >}} This will be useful later on to study the convergence with respect to the number of ''k''-points. The other four lines define the shifts, one per line, expressed in reduced coordinates of the Brillouin zone. Alternatively, one can also define the reciprocal-space mesh by explicitly setting the position and weight of each ''k''-point using the {{< variable "KPoints" >}} or {{< variable "KPointsReduced" >}} variables.

* {{< code-inline >}}{{< variable "KPointsUseSymmetries" >}} = yes{{< /code-inline >}}: this variable controls if symmetries are used or not. When symmetries are used, the code shrinks the Brillouin zone to its irreducible portion and the effective number of ''k''-points is adjusted. 

* {{< code-inline >}}{{< variable "Output" >}} = dos{{< /code-inline >}}: we ask the code to output the density of states.

Here we have taken the value of the grid spacing to be 0.5 bohr. Although we will use this value throughout this tutorial, remember that in a real-life calculation the convergence with respect to the grid spacing must be performed for all quantities of interest.

Note that for periodic systems the default value for the {{< code-inline >}}{{< variable "BoxShape" >}}{{< /code-inline >}} variable is {{< code parallelepiped >}}, although in this case the name can be misleading, as the actual shape also depends on the lattice vectors. This is the only box shape currently available for periodic systems.

## Output
Now run {{< octopus >}} using the above input file. Here are some important things to note from the output.

{{< code-block >}}
#include_file doc/tutorials/periodic_systems/getting_started_with_periodic_systems/1.start/Space.txt
{{< /code-block >}}
This tells us that out system is indeed being treated as periodic in 3 dimensions.

{{< code-block >}}
#include_file doc/tutorials/periodic_systems/getting_started_with_periodic_systems/1.start/Grid.txt
{{< /code-block >}}

Here {{< octopus >}} outputs some information about the cell in real and reciprocal space.

{{< code-block >}}
#include_file doc/tutorials/periodic_systems/getting_started_with_periodic_systems/1.start/Symmetries.txt
{{< /code-block >}}

This block tells us about the space-group and the symmetries found for the specified structure.

{{< code-block >}}
#include_file doc/tutorials/periodic_systems/getting_started_with_periodic_systems/1.start/Lattice.txt
{{< /code-block >}}
Here Octopus outputs some information about the unit cell in real and reciprocal space.

{{< code-block >}}
#include_file doc/tutorials/periodic_systems/getting_started_with_periodic_systems/1.start/k-points.txt
{{< /code-block >}}

Next we get the list of the ''k''-points in reduced coordinates and their weights. Since symmetries are used, only two ''k''-points are generated. If we had not used symmetries, we would have 32 ''k''-points instead.

The rest of the output is much like its non-periodic counterpart. After a few iterations the code should converge:

{{< code-block >}}
#include_file doc/tutorials/periodic_systems/getting_started_with_periodic_systems/1.start/last_iter.txt
#include_file doc/tutorials/periodic_systems/getting_started_with_periodic_systems/1.start/footer.txt
{{< /code-block >}}

As usual, the {{< file "static/info" >}} file contains the most relevant information concerning the calculation. Since we asked the code to output the density of states, we also have a few new files in the {{< file "static" >}} directory:

* {{< file "dos-XXXX.dat" >}}: the band-resolved density of states (DOS);
* {{< file "total-dos.dat" >}}: the total DOS (summed over all bands);
* {{< file "total-dos-efermi.dat" >}}: the Fermi Energy in a format compatible with {{< file "total-dos.dat" >}}.

Of course you can tune the output type and format in the same way you do in a finite-system calculation.

## Convergence in k-points

Similar to the convergence in spacing, a convergence must be performed for the sampling of the Brillouin zone. To do this one must try different numbers of ''k''-points along each axis in the Brillouin zone. This can easily be done by changing the value of the <tt>nk</tt> auxiliary variable in the previous input file. You can obviously do this by hand, but this is something that can also be done with a script. Here is such a script, which we will name {{< file "kpts.sh" >}}:

```bash
#include_file doc/tutorials/periodic_systems/getting_started_with_periodic_systems/2.k-points/kpts.sh
```
 
After running the script ({{< code "source kpts.sh" >}}), you should obtain a file called {{< file "kpts.log" >}} that should look like this:

{{< code-block >}}
#include_file doc/tutorials/periodic_systems/getting_started_with_periodic_systems/2.k-points/kpts.log
{{< /code-block >}}

As you can see, the total energy is converged to within 0.0001 hartree for {{< code "nk = 6" >}}.

You can now play with an extended range, e.g. from 2 to 12. You should then see that the total energy is converged to less than a micro hartree. 

## Band-structure

We now proceed with the calculation of the band-structure of Si. In order to compute a band-structure, we must perform a non-self-consistent calculation, using the density of a previous ground-state calculation. So the first step is to obtain the initial ground-state. To do this, rerun the previous input file, but changing the number of k-points to {{< code "nk=6" >}}. Next, modify the input file such that it looks like this:

{{< code-block >}}
#include_input doc/tutorials/periodic_systems/getting_started_with_periodic_systems/4.band-structure/inp
{{< /code-block >}}

Here are the things we changed:

* {{< code-inline >}}{{< variable "CalculationMode" >}} = unocc{{< /code-inline >}}: we are now performing a non-self-consistent calculation, so we use the {{< code unoccupied >}} calculation mode;

* {{< code-inline >}}{{< variable "ExtraStates" >}} = 10{{< /code-inline >}}: this is the number of unoccupied bands to calculate;

* {{< code-inline >}}{{< variable "ExtraStatesToConverge" >}} = 5{{< /code-inline >}}: the highest unoccupied states are very hard to converge, so we use this variable to specify how many unoccupied states are considered for the stopping criterion of the non-self-consistent run.

* {{< code-inline >}}{{< variable "KPointsPath" >}}{{< /code-inline >}}: this block is used to specify that we want to calculate the band structure along a certain path in the Brillouin zone. This replaces the {{< variable "KPointsGrid" >}} block. The first row describes how many ''k''-points will be used to sample each segment. The next rows are the coordinates of the ''k''-points from which each segment starts and stops. In this particular example, we chose the following path: L-Gamma, Gamma-X, X-Gamma using a sampling of 10-10-15 ''k''-points.

* {{< code-inline >}}{{< variable "KPointsUseSymmetries" >}} = no{{< /code-inline >}}: we have turned off the use of symmetries.

After running {{< octopus >}} with this input file, you should obtain a file named {{< file "bandstructure" >}} inside the {{< file "static" >}} directory. This is how the first few lines of the file should look like:

{{< code-block >}}
#include_file doc/tutorials/periodic_systems/getting_started_with_periodic_systems/4.band-structure/bandstructure-head.txt
...
{{< /code-block >}}


The first column is the coordinate of the ''k''-point along the path. The second, third, and fourth columns are the reduced coordinates of the ''k''-point. The following columns are the eigenvalues for the different bands. In this case there are 14 bands (4 occupied and 10 unoccupied).

Below you can see the plot of the band structure. 
#include_eps doc/tutorials/periodic_systems/getting_started_with_periodic_systems/4.band-structure/Si_bandstructure.eps caption="Band structure of bulk silicon. The zero of energy has been shifted to the maximum of the occupied bands."

This plot shows the occupied bands (purple) and the first 5 unoccupied bands (green). If you are using gnuplot, you can obtain a similar plot with the following command:

```text
plot for [col=5:5+9] 'static/bandstructure' u 1:(column(col)) w l notitle ls 1
```

Note that when using the {{< variable "KPointsPath" >}} input variable, {{< octopus >}} will run in a special mode, and the restart information of the previous ground-state calculation will not be altered in any way. The code informs us about this just before starting the unoccupied states iterations:

```text
Info: The code will run in band structure mode.
     No restart information will be printed.
```

{{< tutorial-footer >}}










