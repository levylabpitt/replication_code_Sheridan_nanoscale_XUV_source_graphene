---
title: "DFT+U"
series: "Tutorials"
author: "Nicolas Tancogne-Dejean"
system_types: "bulk"
description: "how to add the on-site Hubbard interaction U"
---


The objective of this tutorial is to give an idea of how DFT+U in {{< octopus >}} works. As a prototypical example, we will consider bulk NiO in its anti-ferromagnetic configuration.

Details on the implementation of the DFT+U can be found in Ref.[^footnote-1]
.

##  Bulk NiO with PBE  

We will start by setting up an input file for the system we want to simulate, but without DFT+U, as there are a few important things to consider in order to do this correctly. 

###  Input  

The input file we will use is the following one:

{{< code-block >}}
#include_input doc/tutorials/periodic_systems/dft_u/pbe/inp
{{< /code-block >}}


Let's now look more in detail at some of the other input options:

* {{< variable "LatticeVectors" >}} and {{< variable "LatticeParameters" >}}: Here we have decided to neglect the small lattice distortion of bulk NiO in its anti-ferromagnetic configuration and consider only its cubic cell. Note that, as we are interested in the antiferromagnetic order, the primitive cell is doubled along the last lattice vector.

* {{< variable "Species" >}} and {{< variable "ReducedCoordinates" >}}: Octopus uses the {{<code "spglib">}}  library for detecting crystal symmetries, but {{<code "spglib">}} cannot detect the magnetic space group of a system at the present time. Therefore, in order to obtain the correct symmetries, we need to explicitly make the atoms different. To do this, we defined two species, {{<code "Ni1">}} and {{<code "Ni2">}}, which are both Ni atoms, but having different magnetic moments.

* {{<code-inline>}} {{<variable "KPointsUseTimeReversal">}} = no{{</code-inline>}}: Although we obviously want to use the ''k''-point symmetries in our calculation, we need to deactivate the use of time-reversal symmetry, as this is currently not implemented when doing a spin-polarized calculation.

* {{<code-inline>}}{{<variable "GuessMagnetDensity" >}}{{</code-inline>}} and {{<code-inline>}}{{<variable "AtomsMagnetDirection" >}}{{</code-inline>}}: An initial guess is added to break the system symmetries and help the convergence.

We also added a few extra state to the calculation. While this is not needed to obtain the ground state, it is interesting to do so in order to get the value of the electronic bandgap, which we might want to compare to the experiment. Note that depending on the system degeneracy, one might need to add more than one extra state to get an estimate of the bandgap.

###  Output  

After running {{< octopus >}}, we can have a look at the output and notice a few things. The information about the symmetries of the system is printed at the start of the calculation:

{{< code-block >}}
#include_input doc/tutorials/periodic_systems/dft_u/pbe/Symmetries.txt
{{< /code-block >}}

Here we see that the symmetry finder found a different spacegroup for our cubic supercell than the one of cubic NiO where all Ni atoms are considered equivalent (i.e. spacegroup 225).


As the calculation is a spin-polarized calculation, during the SCF cycle {{< octopus >}} outputs the total magnetic moment as well as the local magnetic moments around the atoms, obtained by intergrating the magnetization density on a sphere centered around each atom. For the last SCF iteration, it should look like this:

{{< code-block >}}
#include_input doc/tutorials/periodic_systems/dft_u/pbe/Moments.txt
{{< /code-block >}}

As we can see, we indeed obtained an antiferromagnetic state, where the total magnetization is zero and the local magnetic moments are non-zero and of opposite signs. 

Finally, let's have a look at the values of the direct and indirect bandgaps, which can be found in the {{< file "static/info" >}} file:
{{< code-block >}}
#include_input doc/tutorials/periodic_systems/dft_u/pbe/Gaps.txt
{{< /code-block >}}

{{% notice note %}}
Please, be aware that `ik` is a combined index, describing the k-point and the spin projection. In this case, 'ik=2' still refers to the $\Gamma$ point but for down spins (which are degenerate with the up spins).
{{% /notice %}}

It will be instructive to compare these values with the ones obtained with DFT+U.

##  Bulk NiO with DFT+U  

###  Input  

We will now run the same system, but with a Hubbard U correction. For this we modify our previous input file so that it looks like this:

{{< code-block >}}
#include_input doc/tutorials/periodic_systems/dft_u/dft_u/inp
{{< /code-block >}}


The main differences compared to the previous input file are the following ones:

* {{<code-inline>}}{{< variable "DFTULevel" >}} = dft_u_empirical{{</code-inline>}}: This variable specifies the level of DFT+U used. Here we have chosen to use an empirical correction. The other two available options are {{<code-inline>}}dft_u_none{{</code-inline>}}, which is the default and corresponds to no +U correction, and {{<code-inline>}}dft_u_acbn0{{</code-inline>}}, which corresponds to the ab initio U correction based on the ACBN0 functional[^footnote-2].

* {{<code-inline>}}{{< variable "Species" >}}{{</code-inline>}}: As we have chosen to use an empirical correction, we need to specify the value of the effective Hubbard U and the orbitals this will be aplied to. This is done by adding two options to the {{<code-inline>}}{{< variable "Species" >}}{{</code-inline>}} block: {{<code-inline>}}hubbard_l{{</code-inline>}} specifies the orbitals (l=0 for s orbitals, l=1 for p orbitals, ...) and {{<code-inline>}}hubbard_u{{</code-inline>}} is used to set the value of the effective Hubbard U. Here we add a Hubbard U of 5eV on the 3d orbitals (corresponding to the quantum number l=2).

* {{<code-inline>}}{{< variable "Output" >}}{{</code-inline>}}: We ask the code to output the density matrix of the selected localized subspaces, usually called occupation matrix. This will create a file named {{< file "static/occ_matrices" >}}.

###  Output  

After running {{< octopus >}} using the above input, we can look at the output. Some information specific to DFT+U is printed at the start of the calculation:
{{< code-block >}}
#include_input doc/tutorials/periodic_systems/dft_u/dft_u/DFT+U.txt
{{< /code-block >}}

Among other details, this output confirms that we are indeed using 2 sets of d orbitals, each corresponding to the orbitals centered around a Ni atom.

Lets now look at the magnetic moments:
{{< code-block >}}
#include_input doc/tutorials/periodic_systems/dft_u/dft_u/Moments.txt
{{< /code-block >}}

and bandgaps:

{{< code-block >}}
#include_input doc/tutorials/periodic_systems/dft_u/dft_u/Gaps.txt
{{< /code-block >}}

We see that the magnetic moments on the Ni atoms are larger than at the PBE level and that the bandgap is also found to be larger. Clearly, by adding an Hubbard U correction, we have increased the electron localization in the <tt>d</tt> orbitals, which leads to an increase of the local magnetic moment and to the opening of the electronic bandgap.

##  Double counting term   
The DFT+U method, like DFT+DMFT, needs to specify an approximate treatment of the so called double-counting term, i.e., the part of the electron-electron interaction that is already included in the DFT part and would be double counted if this would not be considered. By default, Octopus uses the fully-localized limit (FLL) double counting. It is also possible to use the around mean-field (AMF) double counting.
This is controlled by the variable {{<variable "DFTUDoubleCounting" >}}.

##  Choice of basis  

In this tutorial, we have seen how to use DFT+U with atomic-center orbitals obtained from the pseudopotential. It is also possible to construct a localization subspace from states defined on the mesh, obtained from a different calculation. For instance, one can perform an DFT calculation to obtain a localized state (core state, flat band, ...). Then is it possible to perform a DFT+U calculation based on these state, using the variable called {{<variable "DFTUBasisFromStates" >}}, together with {{<variable "DFTUBasisStates" >}}.
It is also possible to use oct-wannier90 to obtain the Wannier state and then to use it for the localized subspace. 
This method will be presented in a different tutorial.


{{<tutorial-footer>}}

##  References  
<references/>



[^footnote-1]: {{< article title="Self-consistent DFT+U method for real-space time-dependent density functional theory calculations" authors="Nicolas Tancogne-Dejean, Micael J. T. Oliveira, and Angel Rubio" journal="Phys. Rev. B" year="2017" volume="96" pages="245133" doi="10.1103/PhysRevB.96.245133" >}}

[^footnote-2]: {{< article title="Reformulation of $\mathrm{DFT}+U$ as a Pseudohybrid Hubbard Density Functional for Accelerated Materials Discovery" authors="Agapito, Luis A. and Curtarolo, Stefano and Buongiorno Nardelli, Marco" journal="Phys. Rev. X" volume="5" pages="011006" year="2015" doi="10.1103/PhysRevX.5.011006" >}}

