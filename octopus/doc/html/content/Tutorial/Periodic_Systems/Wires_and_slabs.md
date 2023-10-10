---
title: "Wires and slabs"
tutorials: ["Periodic Systems"]
difficulties: "advanced"
theories: ["DFT"]
calculation_modes: ["Ground State", "Unoccupied"]
system_types: ["Chain", "Slab"]
species_types: "Pseudopotentials"
features: "Band Structure"
description: "Use the flexibility of the real-space grid to treat systems that are periodic in only one or two dimensions."
series: "Tutorials"
Weight: 2
---


In this tutorial we will explain how to use the flexibility of the real-space grid to treat systems that are periodic in only one or two dimensions. As examples we will use a Na chain and a hexagonal boron nitride (h-BN) monolayer. 

## Introduction

In the [Periodic systems](../Periodic systems) tutorial, we saw that the {{< variable "PeriodicDimensions" >}} input variable controls the number of dimensions to be considered as periodic. In that tutorial we only considered the case <tt>{{< variable "PeriodicDimensions" >}} = 3</tt>. Let's now see in detail the different cases:

* <tt>{{< variable "PeriodicDimensions" >}} = 0</tt> (which is the default) gives a finite system calculation, since Dirichlet zero boundary conditions are used at all the borders of the simulation box;

* <tt>{{< variable "PeriodicDimensions" >}} = 1</tt> means that only the ''x'' axis is periodic, while in all the other directions the system is confined. This value must be used to simulate, for instance, a single infinite wire.

* <tt>{{< variable "PeriodicDimensions" >}} = 2</tt> means that both ''x'' and ''y'' axis are periodic, while zero boundary conditions are imposed at the borders crossed by the ''z'' axis. This value must be used to simulate, for instance, a single infinite slab.

* <tt>{{< variable "PeriodicDimensions" >}} = 3</tt> means that the simulation box is a primitive cell for a fully periodic infinite crystal. Periodic boundary conditions are imposed at all borders.


It is important to understand that performing, for instance, a {{< variable "PeriodicDimensions" >}}<tt> = 1</tt> calculation in {{< octopus >}} is not quite the same as performing a {{< variable "PeriodicDimensions" >}}<tt> = 3</tt> calculation with a large supercell. In the infinite-supercell limit the two approaches reach the same ground state, but this does not hold for the excited states of the system.

Another point worth noting is how the Hartree potential is calculated for periodic systems. In fact the discrete Fourier transform that are used internally in a periodic calculation would always result in a 3D periodic lattice of identical replicas of the simulation box, even if only one or two dimensions are periodic. Fortunately {{< octopus >}} includes a clever system to exactly truncate the long-range part of the Coulomb interaction, in such a way that we can effectively suppress the interactions between replicas of the system along non-periodic axes [^cutoff]
. This is done automatically, since the value of {{< variable "PoissonSolver" >}} in a periodic calculation is chosen according to {{< variable "PeriodicDimensions" >}}. See also the variable documentation for {{< variable "PoissonSolver" >}}.

## Atomic positions

An important point when dealing with semi-periodic systems is that the coordinates of the atoms along the aperiodic directions must be centered around 0, as this is the case for isolated systems.
For slabs, this means that the atoms of the slab must be centered around z=0. If this is not the case, the calculation will lead to wrong results.

## Sodium chain

Let us now calculate some bands for a simple single Na chain (i.e. not a crystal of infinite parallel chains, but just a single infinite chain confined in the other two dimensions).

### Ground-state

First we start we the ground-state calculation using the following input file:

{{<code-block>}}
#include_input doc/tutorials/periodic_systems/wires_and_slabs/1.na_chain/1.gs/inp
{{</code-block>}}

Most of these input variables were already introduced in the {{<tutorial "Periodic systems" "Periodic systems">}} tutorial. Just note that the ''k''-points are all along the first dimension, as that is the only periodic dimension.  

The output should be quite familiar, but there are some noteworthy differences.
```text
#include_file doc/tutorials/periodic_systems/wires_and_slabs/1.na_chain/1.gs/Space.txt
```

Here we see that we are indeed running a 3D system that is periodic along one dimension.
```text
#include_file doc/tutorials/periodic_systems/wires_and_slabs/1.na_chain/1.gs/Lattice.txt
```

Although we specified three lattice parameters in the input file corresponding to three lattice vectors, here the code tells us that it's only using one lattice vector. This is because the lattice is only periodic along one direction and therefore it is fully determined by the first vector. The other vectors are ignored when generating the Bravais lattice, but they are still used for two other purposes. First, they allow to specify atomic positions in reduced coordinates, which is useful when copying the coordinates from other codes. Second, they specify the lenght of the paralleliped box used to generate the real-space grid. This is confirmed a few lines later: 
```text
#include_file doc/tutorials/periodic_systems/wires_and_slabs/1.na_chain/1.gs/Grid.txt
```

Note how the lengths along 'y' and 'z' are the same as the lenghts of the corresponding lattice parameters specified in the input file. 

```text
#include_file doc/tutorials/periodic_systems/wires_and_slabs/1.na_chain/1.gs/Hartree.txt
```
Finally, this piece of output confirms that the code is indeed using a cutoff for the calculation of the Hartree potential, as mentioned in the Introduction. The cutoff used is a cylindrical cutoff and, by comparing the FFT grid dimensions with the size of the simulation box, we see that the Poisson solver is using a supercell doubled in size in the y and z directions.

At this point you might want to play around with the number of k-points until you are sure the calculation is converged. 

### Band structure

We now modify the input file in the following way:


{{<code-block>}}
#include_input doc/tutorials/periodic_systems/wires_and_slabs/1.na_chain/2.unocc/inp
{{</code-block>}}
and run the code. You might notice the comments on the the LCAO in the output. What's going on? Why can't a full initialization with LCAO be done?

You can now plot the band structure using the data from the {{< file "static/bandstructure" >}} file, just like in the {{<tutorial "Periodic systems" "Periodic systems">}} tutorial.



In the introduction we mentioned that performing a calculation with <tt>{{< variable "PeriodicDimensions" >}} = 1</tt> is not quite the same as performing a <tt>{{< variable "PeriodicDimensions" >}}= 3</tt> calculation with a large supercell. Lets now check this. Re-run both the ground-state and the unoccupied calculations, but setting <tt>{{< variable "PeriodicDimensions" >}} = 3</tt> in the above input files. Before doing so, make sure you copy the {{< file "static/bandstructure" >}} file to a different place so that it is not overwritten (better yet, run the new calculations in a different folder). You can see the plot of the two band structures on the right. More comments on this in ref. [^cutoff].

#include_eps doc/tutorials/periodic_systems/wires_and_slabs/1.na_chain/5.plot/Na_chain_bands_1D_3D.eps caption="Band structure for a infinite chain of Sodium atoms, calculated for a single chain (purple lines), and a 3D-periodic crystal of chains in a supercell (green lines)."

## h-BN monolayer

Hexagonal boron nitride (h-BN) is an insulator widely studied which has a similar structure to graphene. Here we will describe how to get the band structure of an h-BN monolayer. 

### Ground-state calculation
We will start by calculating the ground-state using the following input file:

{{<code-block>}}
#include_input doc/tutorials/periodic_systems/wires_and_slabs/2.h-bn_monolayer/1.gs/inp
{{</code-block>}}

Most of the file should be self-explanatory, but here is a more detailed explanation for some of the choices:

* <tt>{{< variable "PeriodicDimensions" >}} = 2</tt>: A layer of h-BN is periodic in the ''x''-''y'' directions, but not in the ''z'' direction, so there are two periodic dimensions.

* <tt>{{< variable "LatticeParameters" >}}</tt>: Here we have set the bond length to 1.445 Å. The box size in the ''z'' direction is ''2 L'' with ''L'' large enough to describe a monolayer in the vacuum. Remember that one should always check the convergence of any quantities of interest with the box length value.

If you now run the code, you will notice that the cutoff used for the calculation of the Hartree potential is different than for the Sodium chain, as is to be expected:
```text
#include_input doc/tutorials/periodic_systems/wires_and_slabs/2.h-bn_monolayer/1.gs/Hartree.txt
```


### Band Structure
After the ground-state calculation, we will now calculate the band structure. This is the input for the non-self consistent calculation:
{{<code-block>}}
#include_input doc/tutorials/periodic_systems/wires_and_slabs/2.h-bn_monolayer/2.unocc/inp
{{</code-block>}}

```text
 {{< variable "CalculationMode" >}} = gs
 {{< variable "UnitsOutput" >}} = ev_angstrom
 {{< variable "ExperimentalFeatures" >}} = yes
 
 {{< variable "PeriodicDimensions" >}} = 2
 
 {{< variable "Spacing" >}} = 0.20*angstrom
 
 BNlength = 1.445*angstrom
 a = sqrt(3)*BNlength
 L = 40
 %{{< variable "LatticeParameters" >}}
  a | a | L
 %
 
 %{{< variable "LatticeVectors" >}}
   1.0 | 0.0       | 0.0
  -1/2 | sqrt(3)/2 | 0.0
   0.0 | 0.0       | 1.0
 %
  
 %{{< variable "ReducedCoordinates" >}}
  'B' | 0.0 | 0.0 | 0.0
  'N' | 1/3 | 2/3 | 0.0
 % 
 
 {{< variable "PseudopotentialSet" >}}=hgh_lda
 
 {{< variable "LCAOStart" >}}=lcao_states 
 
 {{< variable "ExtraStatesToConverge" >}}  = 4
 {{< variable "ExtraStates" >}}  = 8
 
 %{{< variable "KPointsPath" >}}
  12  | 7   | 12   - Number of k point to sample each path
  0.0 | 0.0 | 0.0  - Reduced coordinate of the 'Gamma' k point
  1/3 | 1/3 | 0.0  - Reduced coordinate of the 'K' k point
  1/2 | 0.0 | 0.0  - Reduced coordinate of the 'M' k point
  0.0 | 0.0 | 0.0  - Reduced coordinate of the 'Gamma' k point
 %
 {{< variable "KPointsUseSymmetries" >}} = no
```

In this case, we chose the following path to calculate the band structure: Gamma-K, K-M, M-Gamma, with a sampling of 12-7-12 ''k''-points.

Below is the resulting plot of the occupied bands and the first four unoccupied bands from the {{< file "static/bandstructure" >}} file.

#include_eps doc/tutorials/periodic_systems/wires_and_slabs/2.h-bn_monolayer/2.unocc/Tutorial_band_structure_HBN.eps caption="Band structure of a monolayer h-BN."


[^cutoff]: {{< article title="Exact Coulomb cutoff technique for supercell calculations" authors="C. A. Rozzi, D. Varsano, A. Marini, E. K. U. Gross, and A. Rubio" journal="Phys. Rev. B" volume="73" pages="205119" year="2006" doi="10.1103/PhysRevB.73.205119" >}}

{{< tutorial-footer >}}


