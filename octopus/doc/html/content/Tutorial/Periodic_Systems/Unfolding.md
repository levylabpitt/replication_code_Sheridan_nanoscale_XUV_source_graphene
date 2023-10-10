---
title: "Band structure unfolding"
tutorials: "Periodic Systems"
author: "Nicolas Tancogne-Dejean"
difficulties: "basic"
difficulties_weight: 1
theories: "DFT"
calculation_modes: ["Ground state", "Unoccupied"]
system_types: "bulk"
species_types: "Pseudopotentials"
features: ["Band structure", "DOS"]
description: "How to get the bandstructure of a supercell."
weight: 4
---



In this tutorial, we look at how to perform a band-structure unfolding using Octopus. This calculation is done in several steps and each one is described below.

##  Supercell ground-state  

The first thing to do is to compute the ground state of a supercell. In this example we will use bulk silicon. The input file is similar to the one used in the {{<tutorial "Periodic_systems/periodic systems" "Getting started with periodic systems">}} tutorial, but in the present case we will use a supercell composed of 8 atoms. Here is the corresponding input file:

{{< code-block >}}
#include_input doc/tutorials/periodic_systems/band_structure_unfolding/1.supercell_ground-state/inp
{{< /code-block >}}

All these variables should be familiar from other tutorials. Now run {{< octopus >}} using this input file to obtain the ground-state of the supercell.

##  Unfolding setup   

After obtaining the ground state of the supercell, we now define the primitive cell on which we want to unfold our supercell and the specific ''k''-points path that we are interested in. To do this, we take the previous input file and add some new lines:

{{< code-block >}}
#include_input_snippet doc/tutorials/periodic_systems/band_structure_unfolding/2.unfolding_setup/inp unfolding
{{< /code-block >}}


Lets see more in detail the input variables that were added:

* <tt>{{< variable "UnfoldMode" >}} = unfold_setup</tt>: this variable instructs the utility <tt>oct-unfold</tt> in which mode we are running. As a first step, we are running in the <tt>unfold_setup</tt> mode, which generates some files that will be necessary for the next steps.

* <tt>{{< variable "UnfoldLatticeParameters" >}}</tt> specifies the lattice parameters of the primitive cell. This variable is similar to {{< variable "LatticeParameters" >}}.

* <tt>{{< variable "UnfoldLatticeVectors" >}}</tt> specifies the lattice vectors of the primitive cell. This variable is similar to {{< variable "LatticeVectors" >}}.

* <tt>{{< variable "UnfoldKPointsPath" >}}</tt> specifies the ''k''-points path. The coordinates are indicated as reduced coordinates of the primitive lattice. This variable is similar to {{< variable "KPointsPath" >}}.

{{% expand "Expand for the complete input file" %}}
{{< code-block >}}
#include_input doc/tutorials/periodic_systems/band_structure_unfolding/2.unfolding_setup/inp
{{< /code-block >}}
{{% /expand %}}

Now run the <tt>oct-unfold</tt> utility. You will obtain two files. The first one ({{< file "unfold_kpt.dat" >}}) contains the list of ''k''-points of the specified ''k''-point path, but expressed in reduced coordinates of the supercell. These are the ''k''-points for which we need to evaluate the wavefunctions of the supercell. The second file ({{< file "unfold_gvec.dat" >}}) contains the list of reciprocal lattice vectors that relate the ''k''-points of the primitive cell to the one in the supercell.

{{% expand "Expand for 'unfold_kpt.dat'" %}}
{{< code-block >}}
#include_input doc/tutorials/periodic_systems/band_structure_unfolding/2.unfolding_setup/unfold_kpt.dat
{{< /code-block >}}
{{% /expand %}}

{{% expand "Expand for 'unfold_gvec.dat'" %}}
{{< code-block >}}
#include_input doc/tutorials/periodic_systems/band_structure_unfolding/2.unfolding_setup/unfold_gvec.dat
{{< /code-block >}}
{{% /expand %}}

##  Obtaining the wavefunctions  

In order to perform the unfolding, we need the wavefunctions in the supercell at specific k-points. These points are described in the {{< file "unfold_kpt.dat" >}} file that we obtained in the previous step. To get the wavefunctions, we need to run a non self-consistent calculation. Here is the corresponding input file:

{{< code-block >}}
#include_input doc/tutorials/periodic_systems/band_structure_unfolding/3.obtaining_the_wavefunctions/inp
{{< /code-block >}}

In this input file we have changed the calculation mode ({{< variable "CalculationMode" >}} = unocc) and replaced all the ''k''-points related variables (%{{< variable "KPointsGrid" >}}, %{{< variable "KPointsPath" >}}, and %{{< variable "KPoints" >}}) by the line:

{{< code-block >}}
  include unfold_kpt.dat
{{< /code-block >}}

When doing a normal band structure calculation one normally also wants to have some extra states, therefore we have used the {{< variable "ExtraStates" >}} and {{< variable "ExtraStatesToConverge" >}} variables, as explained in the {{<tutorial "periodic_systems/periodic_systems" "Getting started with periodic systems">}} tutorial. In this particular case, we are requesting to converge 4 unoccupied states, which all correspond to the first valence band in the folded primitive cell, as the supercell is 4 times larger than the initial cell.

Now run {{< octopus >}} on this input file. 
{{% notice warning %}}
Note that after you have performed this step, you won't be able to start a time-dependent calculation from this folder. Indeed, the original ground-state wavefunctions will be replaced by the ones from this last calculation, which are incompatible. Therefore it might be a good idea to create a backup of the restart information before performing this step.
{{% /notice %}}

Alternatively, you can continue in a different directory. Instead of copying the restart information to the new folder, you can use the input variable
{{<variable "RestartOptions">}} to specify explicitely from which directory the restart files are to be read.
##  Unfolding  

Now that we have computed the states for the required ''k''-points, we can finally compute the spectral function for each of the ''k''-points of the specified ''k''-point path. This is done by changing unfolding mode to be <tt>{{< variable "UnfoldMode" >}} = unfold_run</tt>, so the final input file should look like this:

{{< code-block >}}
#include_input doc/tutorials/periodic_systems/band_structure_unfolding/4.unfolding/inp
{{< /code-block >}}

This will produce the spectral function for the full path ({{< file "static/ake.dat" >}}) and for each individual points of the path ({{< file "static/ake_XXX.dat" >}}). The content of the {{< file "static/ake.dat" >}} file should look like this:
{{< code-block >}}
#include_input doc/tutorials/periodic_systems/band_structure_unfolding/4.unfolding/ake.txt
...
{{< /code-block >}}


The first column is the coordinate along the ''k''-point path, the second column is the energy eigenvalue and the last column is the spectral function. There are several ways to plot the information contained in this file. One possibility is to plot it as a heat map. For example, if you are using <tt>gnuplot</tt>, you can try the following command:

{{< code-block >}}
 plot 'static/ake.dat' u 1:2:3 w image
{{< /code-block >}}

The unfolded bandstructure of silicon is shown below. How does it compare to the bandstructure compute in the {{<tutorial "Periodic_systems/periodic_systems" "Getting started with periodic systems">}}?

#include_eps doc/tutorials/periodic_systems/band_structure_unfolding/4.unfolding/Si_unfolded.eps caption="Unfolded band structure of bulk silicon."


Note that it is possible to change the energy range and the energy resolution for the unfolded band structure. This is done by specifying the variables {{< variable "UnfoldMinEnergy" >}}, {{< variable "UnfoldMaxEnergy" >}}, and {{< variable "UnfoldEnergyStep" >}}. It is also important to note here that the code needs to read all the unoccupied wavefunctions and therefore might need a large amount of memory, so make sure enough memory is available.

{{< tutorial-footer >}}






