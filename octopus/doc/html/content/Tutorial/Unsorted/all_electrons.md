---
title: "All-electron calculations"
section: "Tutrials"
---


Under certain conditions, it is possible to perform all-electron calculations in the Octopus code.
In this tutorial, we explore how to perform an all-electron calculation for a carbon atom.

Here is the minimal input file needed to perform the calculation
{{< code-block >}}
#include_input doc/tutorials/other/all_electrons/inp
{{< /code-block >}}

We employed here a species type called "species_full_delta". The idea behind this species is to put a delta charge on top of a grid point, as we know the corresponding potential due to this point charge.
The value of valence charge then determines the number of electrons in the simulation.
As usual, the grid spacing and the radius of the box need to be converge.

Due to the delta-type nature of this approach, it suffers from an intrinsic limitation which is that the "atom" needs to be placed on top of a grid point. To lift this constrain, Octopus also implements another species type, called "species_full_gaussian", where one needs to additionally specify a width of the Gaussian associated with the Gaussian charge. 

After running the above input file, one obtains the corresponding eigenvalues and total energy
{{< code-block >}}
#include_file doc/tutorials/other/all_electrons/info.txt
{{< /code-block >}}


These values can directly be compared to the values available on the NIST website for all-electron LDA calculation for the C atom (https://www.nist.gov/pml/atomic-reference-data-electronic-structure-calculations/atomic-reference-data-electronic-7-5). 

While the values are in reasonable agreement, there still present large deviations. This is because we employed here a too large grid spacing, uncapable of capturing the rapid change of the core charge around the nucleus. 
By reducing the grid spacing to a smaller value, one can converge the results of Octopus and recover the all-electron results from the NIST database.

The convergence of the total energy versus the grid spacing is illustrated in Fig. 1.

#include_eps doc/tutorials/other/all_electrons/AllElectronSpacing.eps caption="Fig.1: Total energy convergence versus grid spacing for an all-electron calculation for a carbon atom with a species full delta."


{{<tutorial-footer>}}