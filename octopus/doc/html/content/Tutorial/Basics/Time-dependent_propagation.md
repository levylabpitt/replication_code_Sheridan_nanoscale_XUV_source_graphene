---
title: "Time-dependent propagation"
#tags: ["Basic", "Time-dependent", "Molecule", "Pseudopotentials", "DFT", "Laser"]
series: "Tutorials"
tutorials: ["Octopus Basics"]
theories: "DFT"
calculation_modes: "Time-dependent"
system_types: "Molecule"
species_types: "Pseudopotentials"
difficulties: "basic"
difficulties_weight: 1
features: "Laser"
weight: 7
description: "Time evolution of the density of the methane molecule."
---


OK, it is about time to do TDDFT. Let us perform a calculation of the time evolution of the density of the methane molecule. 

## Ground-state

The reason to do first a ground-state DFT calculation is that this will be the initial state for the real-time calculation.

Run the input file below to obtain a proper ground-state (stored in the {{< code restart >}} directory) as from the 
{{< tutorial "basics:Total Energy Convergence" "total energy convergence tutorial" >}}. (Now we are using a more accurate bond length than before.)

{{< code-block >}}
#include_input doc/tutorials/octopus_basics/time-dependent_run/1.gs/inp
{{< /code-block >}}

## Time-dependent

#### Input

Now that we have a starting point for our time-dependent calculation, we modify the input file so that it looks like this:

{{< code-block >}}
#include_input doc/tutorials/octopus_basics/time-dependent_run/2.td/inp
{{< /code-block >}}

Here we have switched on a new run mode, {{< code-inline >}}{{< variable "CalculationMode" >}} = td{{< /code-inline >}} instead of {{< code gs >}}, to do a time-dependent run. We have also added three new input variables:

* {{< variable "TDPropagator" >}}: which algorithm will be used to approximate the evolution operator. {{< octopus >}} has a large number of possible propagators that you can use (see Ref.[^footnote-1] for an overview).


* {{< variable "TDMaxSteps" >}}: the number of time-propagation steps that will be performed.
* {{< variable "TDTimeStep" >}}: the length of each time step.

Note that, for convenience, we have previously defined a couple of variables, {{< code Tf >}} and {{< code dt >}}. We have made use of one of the possible propagators, {{< code aetrs>}}. The manual explains about the possible options; in practice this choice is usually good except for very long propagation where the {{< code etrs >}} propagator can be more stable.

#### Output

Now run {{< octopus >}}. This should take a few seconds, depending on the speed of your machine. The most relevant chunks of the standard output are

{{< code-block >}}
#include_file doc/tutorials/octopus_basics/time-dependent_run/2.td/input-params.txt
{{< /code-block >}}
...
{{< code-block >}}
#include_file doc/tutorials/octopus_basics/time-dependent_run/2.td/td-start.txt
...
#include_file doc/tutorials/octopus_basics/time-dependent_run/2.td/td-end.txt
{{< /code-block >}}

It is worthwhile to comment on a few things:

* We have just performed the time-evolution of the system, departing from the ground-state, under the influence of no external perturbation. As a consequence, the electronic system does not evolve. The total energy does not change (this you may already see in the output file, the third column of numbers), nor should any other observable change. However, this kind of run is useful to check that the parameters that define the time evolution are correct.
* As the evolution is performed, the code probes some observables and prints them out. These are placed in some files under the directory {{< file "td.general" >}}, which should show up in the working directory. In this case, only two files show up, the {{< file "td.general/energy" >}}, and the {{< file "td.general/multipoles" >}} files. The {{< file "td.general/multipoles" >}} file contains a large number of columns of data. Each line corresponds to one of the time steps of the evolution (except for the first three lines, that start with a {{< code "-" >}} symbol, which give the meaning of the numbers contained in each column, and their units). A brief overview of the information contained in this file follows:
** The first column is just the iteration number.
** The second column is the time.
** The third column is the dipole moment of the electronic system, along the $x$-direction: $\\langle \\Phi (t) \\vert \\hat{x} \\vert \\Phi(t)\\rangle = 
\\int\\!\\!{\\rm d}^3r\\; x\\,n(r)$. Next are the $y$- and $z$-components of the dipole moment.
* The file {{< file "td.general/energy" >}} contains the different components of the total energy.
* It is possible to restart a time-dependent run. Try that now. Just increase the value of {{< variable "TDMaxSteps" >}} and rerun {{< octopus >}}. If, however, you want to start the evolution from scratch, you should set the variable {{< variable "FromScratch" >}} to {{< code yes >}}.

## The time step

A key parameter is, of course, the time step, {{< variable "TDTimeStep" >}}. Before making long calculations, it is worthwhile spending some time choosing the largest time-step possible, to reduce the number of steps needed. This time-step depends crucially on the system under consideration, the spacing, on the applied perturbation, and on the algorithm chosen to approximate the evolution operator. 

In this example, try to change the time-step and to rerun the time-evolution. Make sure you are using {{< variable "TDMaxSteps" >}} of at least 100, as with shorter runs an instability might not appear yet. You will see that for time-steps larger than {{< code "0.0024" >}} the propagation gets unstable and the total energy of the system is no longer conserved. Very often it diverges rapidly or even becomes NaN.

Also, there is another input variable that we did not set explicitly, relying on its default value, {{< variable "TDExponentialMethod" >}}. Since most propagators rely on algorithms to calculate the action of the exponential of the Hamiltonian, one can specify which algorithm can be used for this purpose.

You may want to learn about the possible options that may be taken by {{< variable "TDExponentialMethod" >}}, and {{< variable "TDPropagator" >}} -- take a look at the manual. You can now try some exercises:

* Fixing the propagator algorithm (for example, to the default value), investigate how the several exponentiation methods work (Taylor, Chebyshev, and Lanczos). This means finding out what maximum time-step one can use without compromising the proper evolution.

* And fixing now the {{< variable "TDExponentialMethod" >}}, one can now play around with the various propagators.

## Laser fields

Now we will add a time-dependent external perturbation (a laser field) to the molecular Hamiltonian. For brevity, we will omit the beginning of the file, as this is left unchanged. The relevant part of the input file, with the modifications and additions is:

{{< code-block >}}
#include_input_snippet doc/tutorials/octopus_basics/time-dependent_run/3.laser/inp laser
{{< /code-block >}}

#include_eps doc/tutorials/octopus_basics/time-dependent_run/3.laser/Tutorial_TD_Laser.eps caption="Laser field used in the tutorial" 


The most important variables here are the {{< variable "TDExternalFields" >}} block and the associated {{< variable "TDFunctions" >}} block. You should carefully read the manual page dedicated to these variables: the particular laser pulse that we have employed is the one whose envelope function is a cosine.

Now you are ready to set up a run with a laser field. Be careful to set a total time of propagation able to accommodate the laser shot, or even more if you want to see what happens afterwards. You may also want to consult the meaning of the variable {{< variable "TDOutput" >}}.

A couple of important comments:
* You may supply several laser pulses: simply add more lines to the {{< variable "TDExternalFields" >}} block and, if needed, to the {{< variable "TDFunctions" >}} block.
* We have added the {{< code laser >}} option to {{< variable "TDOutput" >}} variable, so that the laser field is printed in the file {{< file "td.general/laser" >}}. This is done immediately at the beginning of the run, so you can check that the laser is correct without waiting. 
* You can have an idea of the response to the field by looking at the dipole moment in the {{< file "td.general/multipoles" >}}. What physical observable can be calculated from the response?
* When an external field is present, one may expect that unbound states may be excited, leading to ionization. In the calculations, this is reflected by density reaching the borders of the box. In these cases, the proper way to proceed is to setup absorbing boundaries (variables {{< variable "AbsorbingBoundaries" >}}, {{< variable "ABWidth" >}} and {{< variable "ABHeight" >}}).

[^footnote-1]: {{< article title="Propagators for the time-dependent Kohn–Sham equations" authors="A. Castro, M.A.L. Marques, and A. Rubio" journal="J. Chem. Phys." volume="121" pages="3425-3433" year="2004" doi="10.1063/1.1774980" >}}


{{< tutorial-footer >}}


