---
title: "Optical spectra from Sternheimer"
#tags: ["Beginner", "Molecule", "Pseudopotentials", "DFT", "Optical Absorption", "Electromagnetic Response", "Sternheimer"]
series: "Tutorials"
tutorials: "Optical Response"
theories: "DFT"
calculation_modes: "Electromagnetic response"
difficulties: "beginner"
sdifficulties_weight: 2
ystem_types: "Molecule"
species_types: "Pseudopotentials"
features: ["Optical absorption", "Sternheimer"]
Weight: 4
description: "Calculate optical spectra in the frequency domain from linear response."
---


We have just seen how to calculate optical spectra in the time domain with a finite perturbation, and in a frequency-domain, linear-response matrix formulation with the Casida equation. Now we will try a third approach, which is in the frequency domain and linear response but rather than using a pseudo-eigenvalue equation as in Casida, uses a self-consistent linear equation, the Sternheimer equation. This approach is also known as density-functional perturbation theory. It has superior scaling, is more efficient for dense spectra, and is more applicable to nonlinear response. One disadvantage is that one needs to proceed one frequency point at a time, rather than getting the whole spectrum at once. We will find we can obtain equivalent results with this approach for the optical spectra as for time propagation and Casida, by calculating the polarizability and taking the imaginary part.

## Ground state

Before doing linear response, we need to obtain the ground state of the system, for which we can use the same input file as for [Optical spectra from Casida](../Optical spectra from Casida), but we will use a tighter numerical tolerance, which helps the Sternheimer equation to be solved more rapidly. Unlike for Casida, no unoccupied states are required. If they are present, they won't be used anyway.

{{< code-block >}}
#include_input doc/tutorials/optical_response/optical_spectra_from_sternheimer/1.gs/inp
{{< /code-block >}}

## Linear response

Change the {{< variable "CalculationMode" >}} to {{< code "em_resp" >}} and add
the lines below to the input file.
{{< code-block >}}
#include_input_snippet doc/tutorials/optical_response/optical_spectra_from_sternheimer/2.em_resp/inp changes
{{< /code-block >}}

The frequencies of interest must be specified, and we choose them based on the what we have seen from the Casida spectrum. The block above specifies 5 frequencies spanning the range 0 to 8 eV (below the resonances) and 9 frequencies spanning the range 10 to 12 eV (where there are peaks). If we didn't know where to look, then looking at a coarse frequency grid and then sampling more points in the region that seems to have a peak (including looking for signs of resonances in the real part of the polarizability) would be a reasonable approach. We must add a small imaginary part ($\eta$) to the frequency in order to be able to obtain the imaginary part of the response, and to avoid divergence at resonances. The resonances are broadened into Lorentzians with this width. The larger the $\eta$, the easier the SCF convergence is, but the lower the resolution of the spectrum.

To help in the numerical solution, we turn off the preconditioner (which sometimes causes trouble here), and use a linear solver that is experimental but will give convergence much faster than the default one.


{{% expand "Full input file" %}}
{{< code-block >}}
#include_input doc/tutorials/optical_response/optical_spectra_from_sternheimer/2.em_resp/inp
{{< /code-block >}}
{{% /expand %}}

## Spectrum

This Perl script can extract the needed info out of the files in the {{< file "em_resp" >}} directory:

```bash
#include_file doc/tutorials/optical_response/optical_spectra_from_sternheimer/2.em_resp/extract.pl
```

Our result is

```text
#include_file doc/tutorials/optical_response/optical_spectra_from_sternheimer/2.em_resp/spectrum.txt
```

##  See also  

{{< tutorial "Unsorted/Sternheimer linear response" "Sternheimer linear response" >}}

{{< tutorial-footer >}}
