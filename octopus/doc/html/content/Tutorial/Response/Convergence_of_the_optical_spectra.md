---
title: "Convergence of the optical spectra"
Weight: 2
#tags: ["Beginner", "Time-dependent", "Molecule", "Pseudopotentials", "DFT", "Optical Absorption", "oct-propagation_spectrum"]
series: "Tutorials"
tutorials: "Optical Response"
theories: "DFT"
calculation_modes: "Time-dependent"
system_types: "Molecule"
species_types: "Pseudopotentials"
features: "Optical Absorption"
utilities: "oct-propagation_spectrum"
difficulties: "beginner"
difficulties_weight: 2
description: "Convergence study of the spectrum with respect to the grid parameters."
---


In the previous {{< tutorial "Response/Optical Spectra from time-propagation" "Optical Spectra from TD" >}} we have seen how to obtain the absorption spectrum of a molecule. Here we will perform a convergence study of the spectrum with respect to the grid parameters using again the methane molecule as an example.

## Convergence with the spacing

In the {{< tutorial "Basics/Total energy convergence" "Total energy convergence tutorial" >}} we found out that the total energy of the methane molecule was converged to within 0.1 eV for a spacing of 0.18 Å and a radius of 3.5 Å. We will therefore use these values as a starting point for our convergence study. Like for the total energy, the only way to check if the absorption spectrum is converged with respect to the spacing is by repeating a series of calculations, with identical input files except for the grid spacing. The main difference in this case is that each calculation for a given spacing requires several steps:

- Ground-state calculation,
- Three time-dependent propagations with perturbations along ''x'', ''y'' and ''z'', 
- Calculation of the spectrum using the {{< file "oct-propagation_spectrum" >}} utility.

As we have seen before, the $T_d$ symmetry of methane means that the spectrum is identical for all directions, so in this case we only need to perform the convergence study with a perturbation along the ''x'' direction.

All these calculations can be done by hand, but we would rather use a script for this boring an repetitive task. We will then need the following two files.

#### inp
This file is the same as the one used in the previous {{< tutorial "Response/Optical Spectra from time-propagation" "Optical Spectra from TD" >}}:

{{< code-block >}}
#include_input doc/tutorials/optical_response/convergence_of_the_optical_spectra/1.spacing/inp
{{< /code-block >}}


#### spacing.sh

And this is the bash script:

```bash
#include_file doc/tutorials/optical_response/convergence_of_the_optical_spectra/1.spacing/spacing.sh
```

This script performs the three steps mentioned above for a given list of spacings. You will notice that we are not including smaller spacings than 0.18 Å. This is because we know from previous experience that the optical spectrum is well converged for spacings larger than the one required to converge the total energy.

Now run the script typing {{< command-line "source spacing.sh" >}}. This will take some time, so you might want to have a break. Once the script is finished running, we need to compare the spectra for the different spacings. Using your favorite plotting software, e.g. {{< command-line "gnuplot" >}}, plot the second column of the {{< file "cross_section_vector" >}} files vs the first column. You can see on the right how this plot should look like. To make the plot easier to read, we have restricted the energy range and omitted some of the spacings. From the plot it is clear that a spacing of 0.24 Å is sufficient to have the position of the peaks converged to within 0.1 eV.

#include_eps doc/tutorials/optical_response/convergence_of_the_optical_spectra/1.spacing/Absorption_spectrum_CH4_spacing.eps caption="Convergence with spacing of methane absorption spectrum."



## Convergence with the radius

We will now study how the spectrum changes with the radius. We will proceed as for the spacing, by performing a series of calculations with identical input files except for the radius. For this we will use the two following files. 

#### inp
This file is the same as the prevous one, with the exception of the spacing and the time-step. As we have just seen, a spacing of 0.24 Å is sufficient, which in turns allows us to use a larger time-step.

{{< code-block >}}
#include_input doc/tutorials/optical_response/convergence_of_the_optical_spectra/2.radius/inp
{{< /code-block >}}


#### radius.sh

And this is the bash script:

```bash
#include_file doc/tutorials/optical_response/convergence_of_the_optical_spectra/2.radius/radius.sh
```

This script works just like the previous one, but changing the radius instead of the spacing. Now run the script ( {{< command "source radius.sh" >}}) and wait until it is finished.

Now plot the spectra from the {{< file "cross_section_vector" >}} files. You should get something similar to the plot on the right. We see that in this case the changes are quite dramatic. To get the first peak close to 10 eV converged within 0.1 eV a radius of 6.5 Å is necessary. The converged transition energy of 9.2 eV agrees quite well with the experimental value of 9.6 eV and with the TDDFT value of 9.25 eV obtained by other authors.[^footnote-1]

#include_eps doc/tutorials/optical_response/convergence_of_the_optical_spectra/2.radius/Absorption_spectrum_CH4_radius.eps caption="Convergence with spacing of methane absorption spectrum."


What about the peaks at higher energies? Since we are running in a box with zero boundary conditions, we will also have the box states in the spectrum. The energy of those states will change with the inverse of the size of the box and the corresponding peaks will keep shifting to lower energies. However, as one is usually only interested in bound-bound transitions in the optical or near-ultraviolet regimes, we will stop our convergence study here. As an exercise you might want to converge the next experimentally observed transition. Just keep in mind that calculations with larger boxes might take quite some time to run.

[^footnote-1]: {{< article title="Time-Dependent Density Functional Theory Calculations of Photoabsorption Spectra in the Vacuum Ultraviolet Region" authors="N. N. Matsuzawa, A. Ishitani, D. A. Dixon, and T. Uda" journal="J. Phys. Chem. A" volume="105" pages="4953–4962" year="2001" doi="10.1021/jp003937v" >}}


{{< tutorial-footer >}}
