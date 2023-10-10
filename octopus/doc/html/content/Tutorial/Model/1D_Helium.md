---
title: "1D Helium"
weight: 3
series: "Tutorials"
tutorials: "Model Systems"
difficulties: "beginner"
difficulties_weight: 2
theories: "Independent particles"
calculation_modes: ["Ground state", "Unoccupied"]
system_types: "Model"
species_types: "User-defined species"
features: "Visualization"
description: "The helium atom in one dimension which also has two electrons."
---


The next example will be the helium atom in one dimension which also has two electrons, just as we used for the harmonic oscillator. The main difference is that instead of describing two electrons in one dimension we will describe one electron in two dimensions. The calculation in this case is not a DFT one, but an exact solution of the Schrödinger equation -- not an exact solution of the helium atom, however, since it is a one-dimensional model.

### Equivalence between two 1D electrons and one 2D electron

To show that we can treat two electrons in one dimension as one electron in two dimensions, lets start by calling $x\,$ and $y\,$ the coordinates of the two electrons. The Hamiltonian would be (note that the usual Coulomb interaction between particles is usually substituted, in 1D models, by the *soft Coulomb* potential, 
$u(x)=(1+x^2)^{(-1/2)}\\,$):

$$
  \\hat{H} = -\\frac{1}{2}\\frac{\\partial^2}{\\partial x^2}
            -\\frac{1}{2}\\frac{\\partial^2}{\\partial y^2}
  +\\frac{-2}{\\sqrt{1+x^2}}+\\frac{-2}{\\sqrt{1+y^2}}+\\frac{1}{\\sqrt{1+(x-y)^2}}.
$$

Instead of describing two electrons in one dimension, however, we may very well think of one electron in two dimensions,
subject to a external potential with precisely the shape given by:

$$
  \\frac{-2}{\\sqrt{1+x^2}}+\\frac{-2}{\\sqrt{1+y^2}}+\\frac{1}{\\sqrt{1+(x-y)^2}}
  \\,,
$$

Since the Hamiltonian is identical, we will get the same result. Whether we regard $x\,$ and $y\,$ as the coordinates of two different particles in one dimension or as the coordinates of the same particle along the two axes in two dimensions is entirely up to us. (This idea can actually be generalized to treat two 2D particles via a 4D simulation in {{< octopus >}} too!) Since it is usually easier to treat only one particle, we will solve the one-dimensional helium atom in two dimensions. We will also therefore get a "two-dimensional wave-function". In order to plot this wave-function we specify an output plane instead of an axis.

### Input

With the different potential and one more dimension the new input file looks like the following

{{< code-block >}}
#include_input doc/tutorials/model_systems/1D_Helium/1.gs/inp
{{< /code-block >}}

For more information on how to write a potential formula expression, see {{< manual "Basics:Input file" "Input file" >}}.

We named the species "helium" instead of "He" because "He" is already the name of a pseudopotential for the actual 3D helium atom.

### Running

The calculation should converge within 14 iterations. As usual, the results are summarized in the {{< file "static/info" >}} file, where you can find

{{< code-block >}}
...
#include_file doc/tutorials/model_systems/1D_Helium/1.gs/info.txt
{{< /code-block >}}

As we are running with non-interacting electrons, the Hartree, exchange and correlation components of the energy are zero. Also the ion-ion term is zero, as we only have one "ion".

###  Unoccupied States  

Now you can do just the same thing we did for the {{< tutorial "Model:1D_Harmonic_Oscillator" "harmonic oscillator">}} and change the unoccupied calculation mode:

{{< code-block >}}
#include_input doc/tutorials/model_systems/1D_Helium/2.unocc/inp
{{< /code-block >}}


We have added extra states and also restricted the wavefunctions to plot ({{< code-inline >}}{{< variable "OutputWfsNumber" >}} = "1-4,6"{{< /code-inline >}}).

The results of this calculation can be found in the file {{< file "static/eigenvalues" >}}. In this case it looks like

{{< code-block >}}
Some of the states are not fully converged!
Criterion =      0.100000E-05

Eigenvalues [H]
 #st  Spin   Eigenvalue      Occupation     Error
   1   --    -2.238257       1.000000      (2.9E-06)
   2   --    -1.815718       0.000000      (7.9E-07)
   3   --    -1.701549       0.000000      (9.7E-07)
   4   --    -1.629240       0.000000      (9.6E-07)
   5   --    -1.608656       0.000000      (9.3E-07)
   6   --    -1.509599       0.000000      (4.1E-07)
{{< /code-block >}}

Apart from the eigenvalues and occupation numbers we asked {{< octopus >}} to output the wave-functions. To plot them, we will use gnuplot. You start it and type

```bash
#include_file doc/tutorials/model_systems/1D_Helium/2.unocc/plot.gp
```

We plot the ground-state, 1st and 2nd excited-state wave-functions. (If you get this, ignore it: {{< code-inline >}} warning: Cannot contour non grid data. Please use "set dgrid3d".{{< /code-inline >}}) Which correspond to singlet and which to triplet states?

#include_eps doc/tutorials/model_systems/1D_Helium/2.unocc/Wfs-st0001.eps caption="Ground-state of He in 1D" 
#include_eps doc/tutorials/model_systems/1D_Helium/2.unocc/Wfs-st0002.eps caption="1st excited-state of He in 1D" 
#include_eps doc/tutorials/model_systems/1D_Helium/2.unocc/Wfs-st0003.eps caption="2nd excited-state of He in 1D" 


### Exercises

* Calculate the helium atom in 1D, assuming that the 2 electrons of helium do not interact (using {{< variable "Dimensions" >}} = 1). Can you justify the differences?

* See how the results change when you change the interaction. Often one models the Coulomb interaction by $1/\sqrt{a^2+r^2}\,$, and fits the parameter $a\,$ to reproduce some required property.

{{< tutorial-footer >}}
