---
title: "Optical spectra of solids"
series: "Tutorials"
tutorials: ["Periodic Systems"]
difficulties: "beginner"
system_types: "bulk"
calculation_modes: "Time-dependent"
description: "How to obtain an optical spectrum for a solid using a time-dependent calculation."
weight: 3
---


In this tutorial we will explore how to compute the optical properties of bulk systems from real-time TDDFT using {{< octopus >}}. Two methods will be presented, one for computing the optical conductivity from the electronic current and another one computing directly the dielectric function using the so-called "gauge-field" kick method. Note that some of the calculations in this tutorial are quite heavy computationally, so they might require to run in parallel over a several cores to finish in a reasonable amount of time.

##  Ground state   

###  Input  

As always, we will start with the input file for the ground state, here of bulk silicon. In this case we will use a primitive cell of Si, composed of two atoms, as used in the [Getting started with periodic systems](../Getting started with periodic systems) tutorial:

{{<code-block >}}
#include_input doc/tutorials/periodic_systems/optical_spectra_of_solids_conductivity/1.gs/inp
{{</code-block >}}

The description of most variables is given in the {{<tutorial "Periodic_systems/periodic_systems" "Getting started with periodic systems">}}. Here we only added one variable:

* {{< variable "SymmetryBreakDir" >}}: this input variable is used here to specify that we want to remove from the list of possible symmetries to be used the ones that are broken by a perturbation along the x axis. Indeed, we will later apply a perturbation along this direction, which will break the symmetries of the system.

###  Output   

Now run {{< octopus >}} using the above input file. The important thing to note from the output is the effect of the {{< variable "SymmetryBreakDir" >}} input option:

```text
#include_file doc/tutorials/periodic_systems/optical_spectra_of_solids_conductivity/1.gs/Symmetries.txt
```


Here {{< octopus >}} tells us that only 4 symmetries can be used, whereas without {{< variable "SymmetryBreakDir" >}}, 24 symmetries could have been used. Since only some symmetries are used, twelve ''k''-points are generated, instead of two ''k''-points if we would have used the full list of symmetries. If we had not used symmetries at all, we would have 32 ''k''-points instead.

##  Optical conductivity  

We now turn our attention to the calculation of the optical conductivity. For this, we will apply to our system a time-dependent vector potential perturbation in the form of a Heaviside step function. This corresponds to applying a delta-kick electric field. However, because electric fields (length gauge) are not compatible with the periodicity of the system, we work here with vector potentials (velocity gauge).

The corresponding input file is given below:
{{<code-block >}}
#include_input doc/tutorials/periodic_systems/optical_spectra_of_solids_conductivity/2.optical_conductivity/inp
{{</code-block >}}


Most variables have already been discussed in the different prior tutorials. Here are the most important considerations:
* {{< variable "TDOutput" >}}: We ask for outputting the total electronic current. Without this, we will not be able to compute the optical conductivity.
* {{< variable "TDExternalFields" >}}: We applied a vector potential perturbation along the ''x'' axis using.
* {{< variable "TDFunctions" >}}: We defined a step function, an used a strength of the perturbation of 0.01. Note that the strength of the perturbation should be small enough, to guaranty that we remain in the linear response regime. Otherwise, nonlinear effects are coupling different frequencies and a an approach based of the "kick" method is not valid anymore.

Note that if you would have forgotten above to specify {{< variable "SymmetryBreakDir" >}}, the code could would have stopped and let you know that your laser field is breaking some of the symmetries of the system:
```text
******************************* FATAL ERROR ***************************
    ** Fatal Error (description follows)
  * In namespace Lasers:
  *--------------------------------------------------------------------
  * The lasers break (at least) one of the symmetries used to reduce the k-points  .
  * Set SymmetryBreakDir accordingly to your laser fields.
***********************************************************************
```

Now run {{< octopus >}} using the above input file. 

After the calculation is finished, you can find in the {{< file "td.general" >}} folder the {{< file "total_current" >}} file, that contains the computed electronic current. From it, we want to compute the optical conductivity. This is done by running the utility {{< file "oct-conductivity" >}}. 
Before doing so, we need to specify the type of broadening that we want to apply to our optical spectrum. In real time, this correspond to applying a damping function to the time signal. For instance, the Lorentzian broadening that we will employ in this tutorial is given by an exponential damping.
To specify it, simply add the line

{{< code-block >}}
  {{< variable "PropagationSpectrumDampMode" >}} = exponential
{{< /code-block >}}


to your input file before running the {{< file "oct-conductivity" >}} utility.
Running the utility produces the file {{< file "td.general/conductivity" >}}, which contains the optical conductivity $\sigma(\omega)$ of our system.
However, note that at the moment the obtained conductivity is not scaled properly, it is normalized to the volume and proportional to the strength of the perturbation.
Hence, in order to obtain the absorption spectrum, we need to compute

$ \Im [ \epsilon(\omega) ] = \frac{4\pi V}{E_0\omega} \Re [ \sigma(\omega)  ] $

where $V$ is the volume of the system (in this case 10.18$^{3}/4$) and $E_0$ is the strength of the kick (here 0.01). Since our perturbation was along the ''x'' direction, we want to plot the ''x'' component of the conductivity. Therefore, in the end, we need to plot the second column of {{< file "td.general/conductivity" >}} multiplied by appropriate factor versus the energy (first column). You can see how the spectrum should look like in Fig. 1. Because we are computing a current-current response, the spectrum exhibits a divergence at low frequencies. This numerical artifact vanishes when more ''k''-points are included in the simulation.

Looking at the optical spectrum, we observe three main peaks. These correspond to the $E_0$, $E_1$, and $E_1'$ critical points of the bandstructure of silicon. Of course, the shape of the spectrum is not converged because we used too few ''k''-points to sample the Brillouin zone. Remember that you should also perform a convergence study of the spectrum with respect with the relevant parameters. In this case those would be the spacing and the ''k''-points.

#include_eps doc/tutorials/periodic_systems/optical_spectra_of_solids_conductivity/3.spectra/Si_abs_conduc.eps caption="Fig. 1: Absorption spectrum of Si obtained from the computed optical conductivity."


##  Dielectric function using the gauge field kick  

Let us now look at a different way to compute the dielectric function of a solid.
The dielectric function can be defined from the link between the macroscopic total electric field (or vector potential) and the external electric field.
Thus, one can solve the Maxwell equation for computing the macroscopic induced field in order to get the total vector potential acting on the system. 

This is the reasoning behind the gauge-field method, as proposed in Ref.[^footnote-1]

The corresponding input file is given below:

{{<code-block >}}
#include_input doc/tutorials/periodic_systems/optical_spectra_of_solids_gauge_field/2.gauge_field/inp
{{</code-block >}}

Note that before running this input file, you need to redo a GS calculation setting ''nk=5'', to have the same ''k''-point grid as in the above input file.

Compared to the previous TD calculation, we have made few changes:
* {{< variable "GaugeVectorField" >}}: This activates the calculation of the gauge field and its output. This replaces the laser field used in the previous method. 
* {{< variable "KPointsGrid" >}}: We used here a much denser ''k''-point grid. This is needed because of the numerical instability of the gauge-field method for very low number of ''k''-points.

After running this input file with {{< octopus >}}, we can run the utility {{< file "oct-dielectric-function" >}}, which produces the output {{< file "td.general/dielectric_function" >}}. In this case we are interested in the imaginary part of the dielectric function along the ''x'' direction, which correspondents to the third column of the file. You can see how the spectra looks like in Fig. 2.
As we can see, below the bandgap, we observe some numerical noise. This is again due to the fact that the spectrum is not properly converged with respect to the number of ''k''-points and vanishes if sufficient ''k''-points are used.
Apart from this, we now obtained a much more converged spectrum, with the main three features starting to resemble the fully-converged spectrum.

#include_eps doc/tutorials/periodic_systems/optical_spectra_of_solids_gauge_field/3.spectra/Eps2_Si.eps caption="Fig. 2. Absorption spectrum of Si obtained using the gauge-field method."

</br></br>

{{< tutorial-footer >}}
##  References   

[^footnote-1]: {{< article title="Real-space, real-time method for the dielectric function" authors="Bertsch, G. F. and Iwata, J.-I. and Rubio, Angel and Yabana, K." journal="Phys. Rev. B" volume="62" pages="7998" year="2000" doi="10.1103/PhysRevB.62.7998" >}}

