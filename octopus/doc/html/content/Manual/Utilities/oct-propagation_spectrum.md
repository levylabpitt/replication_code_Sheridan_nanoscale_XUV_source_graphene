---
title: "oct-propagation spectrum"
#series: "Manual"
---


### NAME 
oct-propagate_spectrum - Calculates the absorption cross section tensor from the results of a time-propagation run.

### SYNOPSIS 
{{< command "oct-propagate_spectrum" >}}

[oct-propagate_spectrum does not read the standard input: all standard input
will be simply ignored. An input file named {{< file "inp" >}} must be present in the
running directory. Also, oct-propagate_spectrum accepts no command line
arguments, since there is not a standard way to do this with Fortran
90.]

### DESCRIPTION 
This program is one of the {{< octopus >}} utilities.

This utility generates the dipole strength function of the given system. Its main input is the td.general/multipoles file. Output is written to a file called spectrum. This file is made of two columns: energy (in eV or a.u., depending on the units specified in the input file), and dipole strength function (in 1/eV, or 1/a.u., idem).

In the input file, the user may set the {{< variable "PropagationSpectrumTransform" >}} (this should be set to “sine” for proper use), the {{< variable "PropagationSpectrumDampMode" >}} (recommended value is “polynomial”, which ensures fulfilling of the N-sum rule), the {{< variable "PropagationSpectrumStartTime" >}}, the {{< variable "PropagationSpectrumEndTime" >}}, the {{< variable "PropagationSpectrumEnergyStep" >}}, and the {{< variable "PropagationSpectrumMaxEnergy" >}}.

{{< manual-foot prev="Manual:External utilities:oct-photoelectron_spectrum" next="Manual:External utilities:oct-run_periodic_table" >}}
---------------------------------------------
