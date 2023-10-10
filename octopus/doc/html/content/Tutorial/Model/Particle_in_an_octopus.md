---
title: "Particle in an octopus"
#tags: ["Beginner", "Ground State", "Model", "User-defined Species", "Independent Particles", "Visualization"]
weight: 4
series: "Tutorials"
tutorials: "Model Systems"
hidden: "yes"
difficulties: "beginner"
difficulties_weight: 2
theories: "Independent particles"
calculation_modes: "Ground state"
species_types: "User-defined species"
features: "Visualization"
description: "2D systems where the shape of the simulation box is defined by an image file."
---


{{< octopus >}} can actually run 2D systems where the shape of the simulation box is defined by what is white in an image file. Here is an example of a "particle in an octopus", in which we have a constant potential and an octopus-shaped quantum dot. To run it, you will need to have built the code with the [optional library GDLIB](https://libgd.github.io).

## Input

For this example we will need two files:

#### {{< file "inp" >}}

{{< code-block >}}
#include_input doc/tutorials/model_systems/particle_in_a_octopus/inp 
{{< /code-block >}}


#### {{< file "Gdlib.png" >}} 

Make this file available in the run directory. You can download it by clicking on the image bellow. It is also available in the {{< file "PREFIX/share/octopus" >}} directory from your {{< octopus >}} installation.

{{< figure src="/images/Gdlib.png" >}}

## Plotting

The wavefunction are obtained by using:

{{< code-block >}}
 %{{< variable "Output" >}}
   wfs
 %
 {{< variable "OutputFormat" >}} = plane_z
{{< /code-block >}}

View it in {{< code "gnuplot" >}} with

```bash
 plot 'static/wf-st0001.z=0' u 1:2:3 linetype palette
```

or

```bash
 splot 'static/wf-st0001.z=0' u 1:2:(0):($3*500) with pm3d
```

Where does the wavefunction localize, and why?

## Exercises
* See how the total energy scales with the size of the system (controlled by the {{< code ff >}} parameter in the input file). How does it compare to the formula for a particle in a box?
* Look at the wavefunctions of the unoccupied states.
* Think of a serious application that would use the ability to define the simulation box by an image!

{{< tutorial-footer >}}
