---
title: "Centering a geometry"
# tags: ["Basic", "Molecule", "oct-center-geom"]
series: "Tutorials"
tutorials: ["Octopus Basics"]
difficulties: "basic"
difficulties_weight: 1
system_types: "Molecule"
utilities: "oct-center-geom"
description: "Translate the center of mass to the origin."
weight: 5
---


Before running an {{< octopus >}} calculation of a molecule, it is always a good idea to run the {{< manual "Utilities:oct-center-geom"  "oct-center-geom" >}} utility, which will translate the center of mass to the origin, and align the molecule, by default so its main axis is along the ''x''-axis. Doing this is often helpful for visualization purposes, and making clear the symmetry of the system, and also it will help to construct the simulation box efficiently in the code. The current implementation in {{< octopus >}} constructs a parallelepiped containing the simulation box and the origin, and it will be much larger than necessary if the system is not centered and aligned. For periodic systems, these considerations are not relevant (at least in the periodic directions).

## Input
For this example we need two files.

#### {{< file "inp" >}}
We need only a very simple input file, specifying the coordinates file.

{{< code-block >}}
#include_input doc/tutorials/octopus_basics/centering_a_geometry/inp
{{< /code-block >}}

#### {{< file "tAB.xyz" >}}

We will use this coordinates file, for the molecule [''trans''-azobenzene](https://en.wikipedia.org/wiki/Azobenzene), which has the interesting property of being able to switch between ''trans'' and ''cis'' isomers by absorption of light.

{{< code-block >}}
#include_file doc/tutorials/octopus_basics/centering_a_geometry/tAB.xyz
{{< /code-block >}}

## Centering the geometry

{{< manual "Utilities:oct-center-geom" "oct-center-geom" >}} is a serial utility, and it will be found in the {{< code "bin" >}} directory after installation of {{< octopus >}}.

When you run the utility, you should obtain a new coordinates file {{< code "adjusted.xyz" >}} for use in your calculations with {{< octopus >}}. The output is not especially interesting, except for perhaps the symmetries and moment of inertia tensor, which come at the end of the calculation:
 
{{< code-block >}}
#include_file doc/tutorials/octopus_basics/centering_a_geometry/Symmetries.txt
{{< /code-block >}}

You can now visualize the original and new coordinates files with {{< code "xcrysden" >}} or your favorite visualization program (''e.g.'' Jmol, Avogadro, VMD, Vesta, etc.), and see what transformations the utility performed. You can see more options about how to align the system at {{< variable "AxisType" >}} and {{< variable "MainAxis" >}}.

{{< tutorial-footer >}}

