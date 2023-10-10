---
title: "Casida Linear Response"
#series: "Manual"
weight: 3
description: "Alternative to linear response"
---


Mark Casida's formulation of Linear-Response TDDFT allows calculations of the excitation energies of a finite system. For small molecules, this is normally the fastest way to calculate them.

To perform a Casida calculation you first need a {{< manual "Calculations/Ground State" "ground-state" >}} calculation and a calculation of unoccupied states; then you run the code with {{< variable "CalculationMode" >}}{{< code "=casida" >}}.

Here are the input files for a linear-response calculation of the nitrogen dimer, found at {{< file "share/testsuite/linear_response/01-casida.*" >}}.

#### Ground state calculation:

The first step is to converge the ground-state:

{{< expand "Ground State input file" >}}
{{< code-block >}}
#include_input testsuite/linear_response/01-casida.01-gs.inp
{{< /code-block >}}
{{< /expand >}}

#### Unoccupied states

Once the ground state calculation is converged, it is necessary to add and converge more unoccupied states.

{{< expand "Unoccupied States input file" >}}
{{< code-block>}}
#include_input testsuite/linear_response/01-casida.0?-unocc.inp
{{< /code-block >}}
{{< /expand >}}

#### Casida calculation for excited states

{{< expand "Casida input file" >}}
{{< code-block>}}
#include_input testsuite/linear_response/01-casida.0?-casida.inp
{{< /code-block >}}
{{< /expand >}}

See also the {{< tutorial "Response/Optical_spectra_from_Casida" "tutorial on Casida calculations">}}.

{{< manual-foot prev="Calculations:Time-Dependent" next="Calculations:Linear_Response" >}}
