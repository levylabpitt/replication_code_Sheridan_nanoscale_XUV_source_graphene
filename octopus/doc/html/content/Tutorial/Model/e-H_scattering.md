---
title: "e-H scattering"
weight: 5
series: "Tutorials"
tutorials: "Model Systems"
difficulties: "advanced"
difficulties_weight: 2
theories: "DFT"
calculation_modes: "Time dependent"
species_types: "User-defined species"
features: "Visualization"
description: "Electron wave packet scattering on an hydrogen atom"
---


In this tutorial, we will show how to model the problem of an electron wavepacket scattering on an hydrogen atom using {{< octopus >}}. In order to speed up the calculations and make it easier to plot the different quantities of interest, we will do the calculation in 1D. It should be straightforward to change this example to simulate a "real" 3D hydrogen atom.

##  Ground state  

The first thing to do is to tell {{< octopus >}} what we want it to do. We will start with an input file containing only the description of a 1D hydrogen atom, modeled by a soft Coulomb potential:

{{< code-block >}}
#include_input doc/tutorials/model_systems/e-H_scattering/1.gs/inp
{{< /code-block >}}

These input variables should already be familiar as they have been explained in previous tutorials.
Here we simply note the following:
* {{< variable "Spacing" >}}: We employed a large box. This is needed for the time-dependent simulation, and, as we will see, is in fact not large enough for correctly describing our e-H problem for even 1 fs.
* {{< variable "ExtraStates" >}}: We requested one extra state. This is not really needed for the ground-state calculation, but is crucial for the time-dependent run bellow, as there we will replace this unoccupied state by the electron wavepacket.
* {{< variable "Occupations" >}}: We have fixed the occupations to be one in the first spin channel. This is needed for the TD run, as we will see, and also because otherwise the code might get to a different ground state than the one we are interested.

##  Output  
Now one can execute this file by running {{< octopus >}}. Here are a few things worthy of a closer look in the standard output.
We are doing an LDA calculation, as we want to investigate the role of electron-electron interaction here, and, as a possible application, benchmark the performance of exchange-correlation functionals.
As we are doing a one dimensional calculation, we see that {{< octopus >}} has selected the corresponding 1D version of the LDA functional:


{{< code-block >}}
#include_file doc/tutorials/model_systems/e-H_scattering/1.gs/theory_level.txt
{{< /code-block >}}

##  e-H scattering  

We are now ready to perform a simulation of the e-H scattering problem.
For this, we use the following input file:

{{< code-block >}}
#include_input doc/tutorials/model_systems/e-H_scattering/2.td/inp
{{< /code-block >}}

In order to introduce the electron wavepacket to scatter on the hydrogen atom, we have done the following change:
* {{< variable "ExcessCharge" >}}: We introduce an extra electron in our simulation by setting this variable to -1.
* {{< variable "Occupations" >}}: We force the extra electron to be on the second state of the first spin channel.
* {{< variable "RestartFixedOccupations" >}}<tt> = no</tt>: Note, that we need to override the default behaviour to keep the occupations from the restart file.
* {{< variable "UserDefinedStates" >}}: We are replacing the originally unoccupied state obtained after the GS calculation by a user defined state, which is a Gaussian of the form $
\phi(x) = e^{-\alpha (x-x_0)^2+i p (x-x_0)}, $ where $\alpha$ is taken to be $\alpha=0.1$ here, $x_0$ is the initial position of the wavepacket, taken to be at $x_0=-10$ a.u., and $p$ its velocity, taken to be $p=-1.5$ a.u..


The result of the calculation is shown in Fig. 1. What we are interested in is the density of the up channel, where both electrons are, as well as the corresponding potential. These quantities are found in the file 'output_iter/td.XXXXXXX/density-sp1.y=0,z=0' and  'output_iter/td.XXXXXXX/vxc-sp1.y=0,z=0'.
We note a few interesting details. First, we observe no significant reflection of the wavepacket density on the hydrogen atom. This is a known deficiency of the adiabatic approximation.
Moreover, towards the end of the simulation, both the density and the exchange-correlation potential show strong oscillations on the left-hand side of the simulation.
This is due to the artificial reflection of the electron wavepacket at the border of the simulation box.
This can be easily fixed by double the size of the box ({{< variable "Radius" >}} = 100), and/or by using absorbing boundaries.

### Gnuplot script 
In order to produce the above movie, one can use the following gnuplot script
```bash
#include_file doc/tutorials/model_systems/e-H_scattering/2.td/plot.gp
```

and convert it to a movie using convert
```bash
#include_file doc/tutorials/model_systems/e-H_scattering/2.td/convert.sh
```

{{< figure src="/images/scattering.gif" width="500px" caption="Fig. 1. e-H scattering computed at the ALDA level." >}}

{{<tutorial-footer>}}






