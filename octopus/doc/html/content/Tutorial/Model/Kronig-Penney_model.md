---
title: "Kronig-Penney Model"
section: "/home/lueders/Downloads/Tutorial Kronig-Penney_Model"
series: "Tutorials"
tutorials: "Model Systems"
weight: 5
theories: "Independent particles"
calculation_modes: "Ground state"
system_types: "Chain"
species_types: "User-defined species"
features: "Band structure"
difficulties: "beginner"
difficulties_weight: 2
description: "Calculate the bandstructure for Kronig-Penney Model."
---


The Kronig-Penney model is a 1D system that demonstrates band gaps, which relate to the allowed energies for electrons in a material. In this tutorial we calculate the bandstructure for Kronig-Penney Model. The Kronig-Penney Model has a periodic potential of

$$
V(x) =
\begin{cases}
      V_0 & -b < x < 0 \cr
      0 & 0 < x < a
\end{cases}
$$

Where b is the width of each barrier, and a is the spacing between them. 

## Input
The following input file will be used for the ground state calculation:

{{< code-block >}}
#include_input doc/tutorials/model_systems/kronig_penney_model/1.gs/inp
{{< /code-block >}}

#include_eps doc/tutorials/model_systems/kronig_penney_model/1.gs/wavefunctions.eps caption="The first two wavefunctions plotted alongside the potential."

## Bandstructure
To calculate the bandstructure simply change the {{< variable "CalculationMode" >}} to unocc.

{{< code-block >}}
#include_input doc/tutorials/model_systems/kronig_penney_model/2.unocc/inp
{{< /code-block >}}


#include_eps doc/tutorials/model_systems/kronig_penney_model/2.unocc/bandstructure.eps caption="The band structure for Kronig-Penney Model."


To plot the bandstructure, we will use the same command from the {{< tutorial "Periodic_Systems" "Periodic Systems" >}} (assuming you are using gnuplot).

```bash
#include_file doc/tutorials/model_systems/kronig_penney_model/2.unocc/plot.gp
```

Reference:

Sidebottom DL. Fundamentals of condensed matter and crystalline physics: an introduction for students of physics and materials science. New York: Cambridge University Press; 2012.



{{< tutorial-footer >}}





