---
title: "Optical spectra from Casida"
#tags: ["Beginner", "Unoccupied", "Casida", "Molecule", "Pseudopotentials", "DFT", "Optical Absorption", "oct-casida_spectrum"]
series: "Tutorials"
tutorials: "Optical Response"
theories: "DFT"
calculation_modes: "Time-dependent"
difficulties: "beginner"
difficulties_weight: 2
system_types: "Molecule"
species_types: "Pseudopotentials"
features: "Optical absorption"
utilities: "oct-casida_spectrum"
Weight: 3
#series: "Tutorial"
description: "Calculate the absorption spectrum of methane using Casida's equations"
---


In this tutorial we will again calculate the absorption spectrum of methane, but this time using Casida's equations.


## Ground-state

Once again our first step will be the calculation of the ground state. We will use the following input file:

{{< code-block >}}
#include_input doc/tutorials/optical_response/optical_spectra_from_casida/1.gs/inp
{{< /code-block >}}

Note that we are using the values for the spacing and radius that were found in the {{< tutorial "Response/Convergence of the optical spectra" "Convergence_of_the_optical_spectra tutorial" >}} to converge the absorption spectrum.

## Unoccupied States

The Casida equation is a (pseudo-)eigenvalue equation written in the basis of particle-hole states. This means that we need both the occupied states -- computed in the ground-state calculation -- as well as the unoccupied states, that we will now obtain, via a non-self-consistent calculation using the density computed in {{< code gs >}}. The input file we will use is

{{< code-block >}}
#include_input doc/tutorials/optical_response/optical_spectra_from_casida/2.unocc/inp
{{< /code-block >}}

Here we have changed the {{< variable "CalculationMode" >}} to {{< code unocc>}} and added 10 extra states by setting the {{< variable "ExtraStates" >}} input variable.

By running {{< octopus >}}, you will obtain the first 10 unoccupied states (do not forget to run the ground-state calculation first). The solution of the unoccupied states is controlled by the variables in section {{< variable "Eigensolver" >}}. You can take a look at the eigenvalues of the unoccupied states in the file {{< file "static/eigenvalues" >}}:

{{< code-block >}}
#include_file doc/tutorials/optical_response/optical_spectra_from_casida/2.unocc/eigenvalues.txt
{{< /code-block >}}

##  Casida calculation  

Now modify the {{< variable "CalculationMode" >}} to {{< code casida>}} and rerun {{< octopus >}}. Note that by default {{< octopus >}} will use all occupied and unoccupied states that it has available. 

Sometimes, it is useful not to use all states. For example, if you have a molecule with 200 atoms and 400 occupied states ranging from -50 to -2 eV, and you are interested in looking at excitations in the visible, you can try to use only the states that are within 10 eV from the Fermi energy. You could select the states to use with {{< variable "CasidaKohnShamStates" >}} or {{< variable "CasidaKSEnergyWindow" >}}.

A new directory will appear named {{< file "casida" >}}, where you can find the file {{< file "casida/casida" >}}:

{{< code-block >}}
#include_file doc/tutorials/optical_response/optical_spectra_from_casida/3.casida/casida.txt
...
{{< /code-block >}}

The \<x\>, \<y\> and \<z\> are the transition dipole moments:

$$
  \<x\> = \<\\Phi\_0|x|\\Phi\_I\>
  \\,;\\qquad
  \<y\> = \<\\Phi\_0|y|\\Phi\_I\>
  \\,;\\qquad
  \<z\> = \<\\Phi\_0|z|\\Phi\_I\>
$$

where $\Phi_0$ is the ground state and $\Phi_I$ is the given excitation. The
oscillator strength is given by:

$$
  f\_I = \\frac{2 m\_e}{3 \\hbar^2} \\omega\_I \\sum\_{n\\in x,y,z} |\<\\Phi\_0|n|\\Phi\_I\>|^2\\,
$$

as the average over the three directions. The optical absorption spectrum can be given as the "strength function",
which is

$$
  S(\\omega) = \\sum\_I f\_I \\delta(\\omega-\\omega\_I)\\,
$$

Note that the excitations are degenerate with degeneracy 3. This could already be expected from the $T_d$ symmetry of methane.

{{% notice note %}}
Note that within the degenerate subspaces, there is some arbitrariness (possibly dependent on the details of your compilation and machine) in the linear combinations of transitions. Therefore, you should not be concerned if you do not have the same results for the components of the transition dipole moments (above) and analysis of the excitations (below). However, the energies and resulting spectra should agree.
{{% /notice %}}

Further information concerning the excitations can be found in the directory {{< file "casida/casida_excitations" >}}. For example, the first excitation at 9.27 eV is analyzed in the file {{< file "casida/casida_excitations/00001" >}}:

{{< code-block >}}
#include_file doc/tutorials/optical_response/optical_spectra_from_casida/3.casida/casida-excitations.txt
...
{{< /code-block >}}

These files contain basically the eigenvector of the Casida equation. The first two columns are respectively the index of the occupied and the index of the unoccupied state, the third is the spin index (always 1 when spin-unpolarized), and the fourth is the coefficient of that state in the Casida eigenvector. This eigenvector is normalized to one, so in this case one can say that 74.7% (0.864<sup>2</sup>) of the excitation is from state 3 to state 5 (one of 3 HOMO orbitals->LUMO) with small contribution from some other transitions.

##  Absorption spectrum  

#include_eps doc/tutorials/optical_response/optical_spectra_from_casida/4.spectrum/Absorption_spectrum_CH4_casida.eps caption="Absorption spectrum of methane calculated with time-propagation and with the Casida equation."

To visualize the spectrum, we need to broaden these delta functions with the utility {{< file "oct-casida_spectrum" >}} (run in your working directory, not in {{< file "casida" >}}). It convolves the delta functions with Lorentzian functions. The operation of {{< file "oct-casida_spectrum" >}} is controlled by the variables {{< variable "CasidaSpectrumBroadening" >}} (the width of this Lorentzian), {{< variable "CasidaSpectrumEnergyStep" >}}, {{< variable "CasidaSpectrumMinEnergy" >}}, and {{< variable "CasidaSpectrumMaxEnergy" >}}. If you run  {{< file "oct-casida_spectrum" >}} you obtain the file {{< file "casida/spectrum.casida" >}}. It contains all columns of {{< file "casida" >}} broadened. If you are interested in the total absorption spectrum, then you should plot the first and fifth columns. You should obtain a picture like the one on the right.

Comparing the spectrum obtained with the time-propagation in the {{< tutorial "Response/Convergence of the optical spectra" "Convergence of the optical spectra tutorial" >}} with the one obtained with the Casida approach using the same grid parameters, we can see that

* The peaks of the time-propagation are broader. This can be solved by either increasing the total propagation time, or by increasing the broadening in the Casida approach.
* The first two peaks are nearly the same. Probably also the third is OK, but the low resolution of the time-propagation does not allow to distinguish the two close peaks that compose it.
* For high energies the spectra differ a lot. The reason is that we only used 10 empty states in the Casida approach. In order to describe better this region of the spectrum we would need more. This is why one should always check the convergence of relevant peaks with respect to the number of empty states.

You probably noticed that the Casida calculation took much less time than the time-propagation. This is clearly true for small or medium-sized systems. However, the implementation of Casida in {{< octopus >}} has a much worse scaling with the size of the system than the time-propagation, so for larger systems the situation may be different. Note also that in Casida one needs a fair amount of unoccupied states which are fairly difficult to obtain.

If you are interested, you may also compare the Casida results against the Kohn-Sham eigenvalue differences in {{< file "casida/spectrum.eps_diff" >}} and the Petersilka approximation to Casida in {{< file "casida/spectrum.petersilka" >}}.

{{< tutorial-footer >}}
