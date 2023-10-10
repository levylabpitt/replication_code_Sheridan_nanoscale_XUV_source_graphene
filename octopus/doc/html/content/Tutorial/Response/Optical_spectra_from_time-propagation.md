---
title: "Optical spectra from time-propagation"
Weight: 1
#tags: ["Beginner", "Time-dependent", "Molecule", "Pseudopotentials", "DFT", "Optical Absorption", "oct-propagation_spectrum"]
difficulties: "beginner"
difficulties_weight: 2
theories: "DFT"
calculation_modes: "Time-dependent"
system_types: "Molecule"
species_types: "Pseudopotentials"
features: "Optical Absorption"
utilities: "oct-propagation_spectrum"
series: "Tutorials"
tutorials: "Optical Response"
description: "Absorption spectrum of a molecule from the explicit solution of the time-dependent Kohn-Sham equations: methane"
---


In this tutorial we will learn how to obtain the absorption spectrum of a molecule from the explicit solution of the time-dependent Kohn-Sham equations. We choose as a test case methane (CH<sub>4</sub>).

## Ground state

Before starting out time-dependent simulations, we need to obtain the ground state of the system. For this we use basically the same {{< file "inp" >}} file as in the 
{{< tutorial "Basics/Total energy convergence" "total energy convergence tutorial" >}}:

{{< code-block >}}
#include_input doc/tutorials/optical_response/optical_spectra_from_time-propagation/1.gs/inp
{{< /code-block >}}

After running {{< octopus >}}, we will have the Kohn-Sham wave-functions of the ground-state in the directory {{< file "restart/gs" >}}. As we are going to propagate these wave-functions, they have to be well converged. It is important not only to converge the energy (that is relatively easy to converge), but the density. 
The default {{< code-inline >}}{{< variable "ConvRelDens" >}} = 1e-06{{< /code-inline >}} is usually enough, though.

## Time-dependent run

To calculate absorption, we excite the system with an infinitesimal electric-field pulse, and then propagate the time-dependent Kohn-Sham equations for a certain time ''T''. The spectrum can then be evaluated from the time-dependent dipole moment.

#### Input

This is how the input file should look for the time propagation. It is similar to the one from the {{< tutorial "Basics/Time-dependent propagation" "Time-dependent propagation tutorial" >}}.

{{< code-block >}}
#include_input doc/tutorials/optical_response/optical_spectra_from_time-propagation/2.td_x/inp_x
{{< /code-block >}}

Besides changing the {{< variable "CalculationMode" >}} to {{< code td >}}, we have added {{< code-inline >}}{{< variable "FromScratch" >}} = yes{{< /code-inline >}}. This will be useful if you decide to run the propagation for other polarization directions (see bellow). For the time-evolution we use again the Approximately Enforced Time-Reversal Symmetry (aetrs) propagator. The time-step is chosen such that the propagation remains numerically stable. You should have learned how to set it up in the tutorial {{< tutorial "Baics/Time-dependent propagation" "time-dependent propagation tutorial" >}}. Finally, we set the number of time steps with the variable {{< variable "TDMaxSteps" >}}. To have a maximum propagation time of 10 $\hbar/{\rm eV}$ we will need around 4350 iterations.

We have also introduced two new input variables to define our perturbation:

* {{< variable "TDDeltaStrength" >}}: this is the strength of the perturbation. This number should be small to keep the response linear, but should be sufficiently large to avoid numerical problems.

* {{< variable "TDPolarizationDirection" >}}: this variable sets the polarization of our perturbation to be on the first axis (''x'').

Note that you will be calculating the singlet dipole spectrum. You can also obtain the triplet by using {{< variable "TDDeltaStrengthMode" >}}, and other multipole responses by using {{< variable "TDKickFunction" >}}. For details on triplet calculations see the {{< tutorial "Response/Triplet excitations" "Triplet excitations tutorial" >}}.

You can now start {{< octopus >}} and go for a quick coffee (this should take a few minutes depending on your machine). Propagations are slow, but the good news is that they scale very well with the size of the system. This means that even if methane is very slow, a molecule with 200 atoms can still be calculated without big problems.

#### Output

The output should be very similar to the one from the {{< tutorial "Basics/Time-dependent propagation" "Time-dependent propagation tutorial" >}}. The main difference is the information about the perturbation:

{{< code-block >}}
#include_input doc/tutorials/optical_response/optical_spectra_from_time-propagation/2.td_x/kick-info.txt
{{< /code-block >}}

You should also get a {{< file "td.general/multipoles" >}} file, which contains the necessary information to calculate the spectrum. The beginning of this file should look like this:

{{< code-block >}}
#include_input doc/tutorials/optical_response/optical_spectra_from_time-propagation/2.td_x/multipoles-head.txt
{{< /code-block >}}

Note how the dipole along the ''x'' direction (forth column) changes in response to the perturbation. 

##  Optical spectra  

In order to obtain the spectrum for a general system one would need to perform three time-dependent runs, each for a perturbation along a different Cartesian direction (''x'', ''y'', and ''z''). In practice this is what you would need to do:
- Set the direction of the perturbation along ''x'' ({{< code-inline >}}{{< variable "TDPolarizationDirection" >}} = 1{{< /code-inline >}}),
- Run the time-propagation,
- Rename the {{< file "td.general/multipoles" >}} file to {{< file "td.general/multipoles.1" >}},
- Repeat the above step for directions ''y'' ({{< code-inline >}}{{< variable "TDPolarizationDirection" >}} = 2{{< /code-inline >}}) and ''z'' ({{< code-inline >}}{{< variable "TDPolarizationDirection" >}} = 3{{< /code-inline >}}) to obtain files {{< file "td.general/multipoles.2" >}} and {{< file "td.general/multipoles.3" >}}.

Nevertheless, the $T_d$ symmetry of methane means that the response is identical for all directions and the absorption spectrum for ''x''-polarization will in fact be equivalent to the spectrum averaged over the three directions. You can perform the calculations for the ''y'' and ''z'' directions if you need to convince yourself that they are indeed equivalent, but you this will not be necessary to complete this tutorial.

Note that {{< octopus >}} can actually use the knowledge of the symmetries of the system when calculating the spectrum. However, this is fairly complicated to understand for the purposes of this introductory tutorial, so it will be covered in the 
{{< tutorial "Response/Use of symmetries in optical spectra from time-propagation" "Use of symmetries in optical spectra from time-propagation tutorial" >}}.

### The {{< file "oct-propagation_spectrum" >}} utility

{{< octopus >}} provides an utility called {{< manual "utilities/oct-propagation_spectrum" "oct-propagation_spectrum" >}} to process the {{< file "multipoles" >}} files and obtain the spectrum. 

#### Input
This utility requires little information from the input file, as most of what it needs is provided in the header of the {{< file "multipoles" >}} files. If you want you can reuse the same input file as for the time-propagation run, but the following input file will also work:

{{< code-block >}}
 {{< variable "UnitsOutput" >}} = eV_angstrom
{{< /code-block >}}

Now run the utility. Don't worry about the warnings generated, we know what we are doing!

#### Output

This is what you should get:

{{< code-block >}}
#include_input doc/tutorials/optical_response/optical_spectra_from_time-propagation/3.propagation_spectra/spectrum-out.txt
{{< /code-block >}}

If you used the full input file, you will see some more parser warnings, which indicate the unused variables.
{{< code-block >}}
#include_input doc/tutorials/optical_response/optical_spectra_from_time-propagation/3.propagation_spectra/parser-warnings.txt
{{< /code-block >}}
Here, these arise as the utility does not need them, but warnings like this also could indicate a mis-typed variable name.


You will notice that the file {{< file "cross_section_vector" >}} is created. If you have all the required information (either a symmetric molecule and one multipole file, or a less symmetric molecule and multiple multipole files), {{< file "oct-propagation_spectrum" >}} will generate the whole polarizability tensor, {{< file "cross_section_tensor" >}}. If not enough multipole files are available, it can only generate the response to the particular polarization you chose in you time-dependent run, {{< file "cross_section_vector" >}}. 

### Cross-section vector

Let us look first at {{< file "cross_section_vector" >}}:

```
#include_file doc/tutorials/optical_response/optical_spectra_from_time-propagation/3.propagation_spectra/cross_section_vector-1.txt
```

The beginning of the file just repeats some information concerning the run.

{{< code-block >}}
#include_file doc/tutorials/optical_response/optical_spectra_from_time-propagation/3.propagation_spectra/cross_section_vector-2.txt
{{< /code-block >}}

Now comes a summary of the variables used to calculate the spectra. Note that you can change all these settings just by adding these variables to the input file before running {{< file "oct-propagation_spectrum" >}}.  Of special importance are perhaps {{< variable "PropagationSpectrumMaxEnergy" >}} that sets the maximum energy that will be calculated, and {{< variable "PropagationSpectrumEnergyStep" >}} that determines how many points your spectrum will contain. To have smoother curves, you should reduce this last variable.

{{< code-block >}}
#include_file doc/tutorials/optical_response/optical_spectra_from_time-propagation/3.propagation_spectra/cross_section_vector-3.txt
{{< /code-block >}}

Now comes some information from the sum rules. The first is just the ''f''-sum rule, which should yield the number of active electrons in the calculations. We have 8 valence electrons (4 from carbon and 4 from hydrogen), but the sum rule gives a number that is much smaller than 8! The reason is that we are just summing our spectrum up to 20 eV (see {{< variable "PropagationSpectrumMaxEnergy" >}}), but for the sum rule to be fulfilled, we should go to infinity. Of course infinity is a bit too large, but increasing 20 eV to a somewhat larger number will improve dramatically the ''f''-sum rule. The second number is the static polarizability also calculated from the sum rule. Again, do not forget to converge this number with {{< variable "PropagationSpectrumMaxEnergy" >}}.

{{< code-block >}}
#include_file doc/tutorials/optical_response/optical_spectra_from_time-propagation/3.propagation_spectra/cross_section_vector-4.txt
...
{{< /code-block >}}

Finally comes the spectrum. The first column is the energy (frequency), the next three columns are a row of the cross-section tensor, and the last one is the strength function for this run.

The dynamic polarizability is related to optical absorption cross-section via $\sigma \left( \omega \right) = \frac{4 \pi \omega}{c} \mathrm{Im}\ \alpha \left( \omega \right) $ in atomic units, or more generally $4 \pi \omega \tilde{\alpha}\ \mathrm{Im}\ \alpha \left( \omega \right) $ (where $\tilde{\alpha}$ is the fine-structure constant) or $\frac{\omega e^2}{\epsilon_0 c} \mathrm{Im}\ \alpha \left( \omega \right) $. The cross-section is related to the strength function by $S \left( \omega \right) = \frac{mc}{2 \pi^2 \hbar^2} \sigma \left( \omega \right)$.

### Cross-section tensor

If all three directions were done, we would have four files: {{< file "cross_section_vector.1" >}}, {{< file "cross_section_vector.2" >}}, {{< file "cross_section_vector.3" >}}, and {{< file "cross_section_tensor" >}}. The latter would be similar to:

{{< code-block >}}
#include_file doc/tutorials/optical_response/optical_spectra_from_time-propagation/3.propagation_spectra/cross_section_tensor-head.txt
...
{{< /code-block >}}

#include_eps doc/tutorials/optical_response/optical_spectra_from_time-propagation/6.propagation_spectra_tensor/Absorption_spectrum_CH4.eps caption="Absorption spectrum of methane"

The columns are now the energy, the average absorption coefficient (Tr $\sigma/3$), the anisotropy, and 
then all the 9 components of the tensor. The third number is the spin component, just 1 here since it is unpolarized. The anisotropy is defined as

$$
 (\\Delta \\sigma)^2 = \\frac{1}{3} \\{ 3{\\rm Tr}(\\sigma^2) - \[{\\rm Tr}(\\sigma)\]^2 \\}
$$

The reason for this definition is that it is identically equal to:

$$
 (\\Delta \\sigma)^2 = (1/3) \[ (\\sigma\_1-\\sigma\_2)^2 + (\\sigma\_1-\\sigma\_3)^2 + (\\sigma\_2-\\sigma\_3)^2 \]\\,
$$

where $\{\sigma_1, \sigma_2, \sigma_3\}\,$ are the eigenvalues of $\sigma\,$. An "isotropic" tensor is characterized by having three equal eigenvalues, which leads to zero anisotropy. The more different that the eigenvalues are, the larger the anisotropy is.

If you now plot the absorption spectrum (column 5 vs 1 in {{< file "cross_section_vector" >}}, but you can use {{< file "cross_section_tensor" >}} for this exercise in case you did propagate in all three directions), you should obtain the plot shown on the right. Of course, you should now try to converge this spectrum with respect to the calculation parameters. In particular:

* Increasing the total propagation time will reduce the width of the peaks. In fact, the width of the peaks in this methods is absolutely artificial, and is inversely proportional to the total propagation time. Do not forget, however, that the area under the peaks has a physical meaning: it is the oscillator strength of the transition.
* A convergence study with respect to the spacing and box size might be necessary. This is covered in the next tutorial.

Some questions to think about:
* What is the equation for the peak width in terms of the propagation time? How does the observed width compare to your expectation?
* What is the highest energy excitation obtainable with the calculation parameters here? Which is the key one controlling the maximum energy?
* Why does the spectrum go below zero? Is this physical? What calculation parameters might change this?


{{< tutorial-footer >}}
