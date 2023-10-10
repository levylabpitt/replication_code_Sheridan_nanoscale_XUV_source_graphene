---
title: "DFT+U+V"
series: "Tutorials"
author: "Nicolas Tancogne-Dejean"
difficulties: "advanced"
system_types: "bulk"
description: "how to add the intersite interaction V"

---


This tutorial aims at explaining how to perform DFT+U+V calculations in {{< octopus >}}. This correspond to adding not only an on-site Hubbard interaction U, but also an intersite interaction V.
As a prototypical example, we will consider bulk silicon.

The DFT+U+V method, as well as its performances and implementation details are discussed in Ref.[^footnote-1]

###  Input  

The input file we will use is the following one:

```
#include_input /doc/tutorials/periodic_systems/dft_u_v/inp
```


Most of the variable are the same as in the introduction tutorial on solids [Getting started with periodic systems](../Getting started with periodic systems), and the DFT+U related variables are explained in the tutorial [DFT+U](../DFT+U).
For convenience, we are doing here a ground-state calculation with a k-point grid and a k-point path. This allows to get the bandstructure as a direct output of a ground-state calculation. However, the unoccupied states might not be converged with this approach. To perform a proper band-structure calculation, see the tutorial [Getting started with periodic systems](../Getting started with periodic systems).

Compared to a more conventional DFT+U calculation, we note few differences here:
* {{<variable "AOLoewdin" >}}{{<code " = yes">}}: We want to perform a Löwdin orthonormalization of the atomic orbitals. This is very important to avoid double counting between different orbitals.
* {{<variable "UseAllAtomicOrbitals" >}} {{<code "= yes">}} and {{<variable "SkipSOrbitals" >}} {{<code "= no" >}}: Here we want to use all the atomic orbitals available in the pseudopotentials. In this case, the *s* orbitals are important to obtain a meaningful orthonormalization. Note that in Si, we indeed want to use all atomic orbitals. But in some strongly correlated materials, one might as well only use one set of atomic orbital using the Species block, as done for NiO in the DFT+U tutorial.
* {{<variable "ACBN0IntersiteInteraction" >}} {{<code "= yes">}}: We are asking for the intersite interaction (+V). The number of neighbors for which we are including the intersite interaction is given by a cutoff in real space, specified by the variable {{<variable "ACBN0IntersiteCutoff" >}}. 
Here we are using 7 Bohr, which gives us the first nearest neighbor only. In practice, this value needs to be converged, or needs to be chosen if one wants to include only certain shells of neighbors.


###  Output  

After running {{< octopus >}} using the above input, we can look at the output. 
In the static/info file, we find the direct and indirect bandgap
{{< code-block >}}
#include_input doc/tutorials/periodic_systems/dft_u/pbe/Gaps.txt

  Direct gap at ik=    6 of  0.1346 H
  Indirect gap between ik=    6 and ik=   24 of  0.1148 H
{{< /code-block >}}

As a comparison, without adding the +U+V, one finds at the LDA level a much smaller bandgap
{{< code-block >}}
#include_input doc/tutorials/periodic_systems/dft_u_v/Gaps.txt

  Direct gap at ik=    1 of  0.1041 H
  Indirect gap between ik=    3 and ik=    9 of  0.0810 H
{{< /code-block >}}


The effect of bandgap opening can be understood more precisely from the bandstructure. As one can see, the conduction bands have been rigidly shifted toward higher energies, as expected from more advance calculations such as hybrid functionals or GW calculations.

#include_eps doc/tutorials/periodic_systems/dft_u_v/Si_bandstructure.eps caption="Band structure of bulk silicon from DFT+U+V. The zero of energy has been shifted to the maximum of the occupied bands."


</pre>

{{< tutorial-footer >}}

##  References  
<references/>

[^footnote-1]: {{< article title="Parameter-free hybridlike functional based on an extended Hubbard model: DFT+U+V" authors="N. Tancogne-Dejean, and A. Rubio" journal="Physical Review B" volume="102" pages="155117" year="2020" doi="10.1103/PhysRevB.102.155117" >}}

