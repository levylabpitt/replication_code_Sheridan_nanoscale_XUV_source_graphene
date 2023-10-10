---
title: "ARPES"
difficulties: "advanced"
theories: "Independent particles"
system_types: "semi-periodic"
species_types: "User defined potential"
description: "This tutorial aims at introducing the user to the basic concepts needed to calculate ARPES with octopus."
series: "Tutorials"
Weight: 20
---

This tutorial aims at introducing the user to the basic concepts needed to calculate ARPES with {{< octopus >}}.
We choose as test case a fictitious non-interacting two-dimensional system with semi-periodic boundary conditions.

## Ground state

Before starting with the time propagation we need to obtain the ground state of the system.

To simulate ARPES simulations with the tSURRF implementation in {{<octopus>}}[^footnote-1] one needs to explicitly model the surface of the material. This requirement often leads to heavy simulations. For this reason, in this tutorial, we illustrate the procedure on a model system. Calculating ARPES with TDDFT on an ab-initio model for a surface involves the same steps. 

### Input

The system we choose is a 2D toy model potential simulating an atomic chain. The dimensionality is selected by telling octopus to use on 2 spatial dimensions (xy) 
with {{<code-inline>}}{{<variable "Dimensions">}} = 2{{</code-inline>}} and by imposing periodic boundary conditions along x with 
{{<code-inline>}}{{<variable "PeriodicDimensions">}} = 1{{</code-inline>}}.


{{< code-block >}}
#include_input doc/tutorials/other/arpes/01-gs/inp
{{< /code-block >}}


Octopus can calculate ARPES on a path in reciprocal space. For this reason we have to specify the path with the {{<variable "KPointsPath">}} 
block already at the ground state level. This is needed to generate the KS wave functions that will be evolved in the time propagation. 
It also provides the bands structure that we will use as reference to cross-check the quality of the spectrum.

Note that, since this is a toy model with independent electrons, {{<code-inline>}}{{<variable "TheoryLevel">}} = independent_particles{{</code-inline>}}, 
we do not need to sample the BZ and the grid in reciprocal space is constituted only by the gamma point, 

{{< code-block >}}
#include_input_snippet doc/tutorials/other/arpes/01-gs/inp kpoints
{{< /code-block >}}



## Time-dependent run

To calculate the ARPES spectrum we excite the system with a laser field to excite electrons into the vacuum, i.e. on the non-periodic dimension. The photoelectron detection probability constituting ARPES is calculated by analyzing the flux of the ionization current trough a surface positioned at a certain distance from the surface of the system. This current is obtained by propagating in time the KS orbitals under the effect of the laser. 

### Input

This is how the input file should look for the time propagation.

{{% expand "full input file" %}}
{{< code-block >}}
#include_input doc/tutorials/other/arpes/02-td/inp
{{< /code-block >}}
{{% /expand %}}


Where we just changed the {{<variable "CalculationMode">}} to {{<code "td">}} to activate the time propagation. 

{{< code-block >}}
#include_input_snippet doc/tutorials/other/arpes/02-td/inp mode
{{< /code-block >}}

We then have to specify the laser field parameters as following.

{{< code-block >}}
#include_input_snippet doc/tutorials/other/arpes/02-td/inp laser_field
{{< /code-block >}}

Since we deal with periodic system we have to specify the light-matter coupling in the so called "velocity gauge" where the field is described by a time-dependent vector potential.
We specify the field as a "vector_potential" polarized along x and with carrier frequency "wpr"

{{< code-block >}}
#include_input_snippet doc/tutorials/other/arpes/02-td/inp vector_potential
{{< /code-block >}}

and a sin^2 envelope function "probe" specified by an analytical expression of time, t.

{{< code-block >}}
#include_input_snippet doc/tutorials/other/arpes/02-td/inp tdfunctions
{{< /code-block >}}

We also specified the total time of the propagation to be exactly synchronized with the pulse switch-off time by setting 
{{<code-inline>}}{{<variable "TDPropagationTime">}} = Tpr{{</code-inline>}}. 
For more details look at the {{<tutorial "basics/time-dependent_propagation" "Time-dependent propagation">}} tutorial.

Photoelectrons ejected from the system by the laser will eventually bounce back from the boundary of the simulation box along the non-periodic dimension (y). We therefore employ absorbing boundary conditions to prevent spurious reflections with the following code block

{{< code-block >}}
#include_input_snippet doc/tutorials/other/arpes/02-td/inp cap
{{< /code-block >}}


Here we employ complex absorbing potential (CAP) boundaries with {{<code-inline>}}{{<variable "AbsorbingBoundaries">}} = cap{{</code-inline>}}, 
specify the CAP parameters {{<variable "ABCapHeight">}} and {{<variable "ABShape">}} to have maximal absorption in the energy region where we expect the photoelectrons.[^footnote-2] 

Finally we specify the parameters for the evaluation of the ARPES spectrum. 

{{< code-block >}}
#include_input_snippet doc/tutorials/other/arpes/02-td/inp arpes
{{< /code-block >}}


Here we tell octopus to use tSURFF with {{<code-inline>}}{{<variable "PhotoElectronSpectrum">}} = pes_flux {{</code-inline >}}, 
specify the position of the analyzing surface at the onset of the absorbing boundaries with 
{{<code-inline>}}{{<variable "PES_Flux_Lsize">}} = Lmin{{</code-inline>}} and define the parameters of the energy grid for the 
final spectrum with the block {{<variable "PES_Flux_EnergyGrid">}}.
The tSURFF method allows to calculate the momentum-resolved photoelectron probability $P({\bf p})$ on an arbitrary grid in momentum space.
The ARPES spectrum is than obtained rewriting $P({\bf p})$ as a function of the total kinetic energy and the momentum parallel to the surface, ${\bf p}\_\parallel$ 
(equivalent to the crystal momentum ${\bf k}$ ) $P({\bf p}\_\parallel={\bf k}, E=({\bf p}\_\parallel + {\bf p}\_\perp)^2/2)$. 
Since a simple cartesian grid in momentum results in a deformed ARPES grid octopus can generate a grid that compensates the deformation such the final ARPES spectrum is in a cartesian grid. 
This grid is generated with {{<code-inline>}}{{<variable "PES_Flux_ARPES_grid">}} = yes{{</code-inline>}}. 
Since this option is true by default for semi-periodic systems we do not have to specify it in the input file.




### Output

If the code runs correctly the standard output should present a Photoelectron section like this:

{{< code-block >}}
*************************** Photoelectron ****************************
Info: Calculating PES using t-surff technique.
Input: [PES_Flux_Shape = pln]
Input: [PES_Flux_Lsize = (8.150,30.00)]
Input: [PES_Flux_Parallelization = pf_none]
Input: [PES_Flux_Momenutum_Grid = cartesian]
Energy grid (Emin, Emax, DE) [H]:  (   1.711,   1.911, 0.10E-01)
Momentum linear grid (Pmin, Pmax, DP) [me*bH(2pi/h)]:  (   1.850,   1.955, 0.14E+00)
Input: [PES_Flux_ARPES_grid = yes]
Number of points with E < p//^2/2 = 0 [of 462]
Input: [PES_Flux_UseSymmetries = no]
Input: [PES_Flux_RuntimeOutput = no]
Info: Total number of surface points = 16
Info: Total number of momentum points = 21
{{< /code-block >}}

## ARPES

The spectral information on the photoelectrons is stored in {{<file "restart/td/pesflux*">}} binary files and can be analyzed in post-processing.

### The {{< command "oct-photoelectron_spectrum" >}} utility

{{<octopus>}} provides an utility called {{<command "oct-photoelectron_spectrum">}} to process the {{<file "pesflux*">}} files and obtain the spectrum.


###  Input  

The utility by default requires no input from the user in calculations designed to obtain an ARPES spectrum. If the td run has been performed with 
{{<code-inline>}}{{<variable "PES_Flux_ARPES_grid">}} = yes{{</code-inline>}} it should have all the information that it needs. 

Just run the utility in the same path where the inp file resides.


### Output

This is what you should get:

{{< code-block >}}
************************** Kpoint selection **************************
Will use a zero-weight path in reciprocal space with the following points
       2    0.000000 |   -0.500000    0.000000 |
       3    0.000000 |   -0.450000    0.000000 |
       4    0.000000 |   -0.400000    0.000000 |
       5    0.000000 |   -0.350000    0.000000 |
       6    0.000000 |   -0.300000    0.000000 |
       7    0.000000 |   -0.250000    0.000000 |
       8    0.000000 |   -0.200000    0.000000 |
       9    0.000000 |   -0.150000    0.000000 |
      10    0.000000 |   -0.100000    0.000000 |
      11    0.000000 |   -0.050000    0.000000 |
      12    0.000000 |    0.000000    0.000000 |
      13    0.000000 |    0.050000    0.000000 |
      14    0.000000 |    0.100000    0.000000 |
      15    0.000000 |    0.150000    0.000000 |
      16    0.000000 |    0.200000    0.000000 |
      17    0.000000 |    0.250000    0.000000 |
      18    0.000000 |    0.300000    0.000000 |
      19    0.000000 |    0.350000    0.000000 |
      20    0.000000 |    0.400000    0.000000 |
      21    0.000000 |    0.450000    0.000000 |
      22    0.000000 |    0.500000    0.000000 |

**********************************************************************

Read PES restart files.
Zenith axis: (      0.00,       0.00,       1.00)
[ 42/ 42] 100%|****************************************|     --:-- ETA

***************** ARPES cut on reciprocal space path *****************
Done
**********************************************************************
{{< /code-block >}}


If all goes well the file is {{<file "PES_ARPES.path">}} will be created containing the ARPES spectrum evaluated on the k-point path.
The spectrum can be visualized with {{<command "gnuplot">}} as a density plot with 

{{< code-block >}}
set pm3d map
sp "PES_ARPES.path" u 1:4:5
{{< /code-block >}}


{{< figure src="images/arpes_2d.png" width="500px" caption="ARPES spectrum for a 2D atomic chain model system." >}}

The ARPES spectrum looks a bit blocky to improve the resolution one has to increase the number of kpoints in "KPointsPath" (requires recomputing the groundstate) and decrease the energy spacing in {{<variable "PES_Flux_EnergyGrid">}}. 
Keep in mind that larger grids implies heavier simulations and longer run times.


Some questions to think about:

* How does the ARPES spectrum compares with the band structure? Plot one on top of the other with {{<command gnuplot>}}.

* How to choose the photoelectron energy range? In the previous example we chose Emin =  wpr - 0.2 and Emax =  wpr in order to see photoelectrons emitted from the valence band. How can we change {{<variable "PES_Flux_EnergyGrid">}} to see the band below? What about conduction band?

* Play around with the laser carrier, "wpr", and pulse envelope, "Tpr". Can you identify the impact of these parameters on the spectrum?



## References

[^footnote-1]: {{< article title="A First-Principles Time-Dependent Density Functional Theory Framework for Spin and Time-Resolved Angular-Resolved Photoelectron Spectroscopy in Periodic Systems" authors="U. De Giovannini, H. Hübener, and A. Rubio" journal="JCTC" volume="13" pages="265" year="2017" doi="10.1021/acs.jctc.6b00897" >}}

[^footnote-2]: {{< article title="Modeling electron dynamics coupled to continuum states in finite volumes with absorbing boundaries" authors="U. De Giovannini, A. H. Larsen, A. Rubio, and A. Rubio" journal="EPJB" volume="88" pages="1" year="2015" doi="10.1140/epjb/e2015-50808-0" >}}


{{< tutorial-footer >}}


