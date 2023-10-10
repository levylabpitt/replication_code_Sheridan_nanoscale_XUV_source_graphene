---
title: "Total energy convergence"
#tags: ["Basic", "Ground State", "Molecule", "Pseudopotentials", "DFT", "Total Energy"]
series: "Tutorials"
tutorials: ["Octopus Basics"]
difficulties: "basic"
difficulties_weight: 1
theories: "DFT"
calculation_modes: "Ground state"
system_types: "Molecule"
species_types: "Pseudopotentials"
features: "Total energy"
description: "Make sure that the results are converged"
weight: 3
---


In this tutorial we will show how to converge a quantity of interest (the total energy) with respect to the parameters that define the real-space grid. For this we will use two examples, the Nitrogen atom from the {{< tutorial "Basics:Basic_input_options" "Basic input options" >}} tutorial, and a methane molecule.

## Nitrogen atom: finding a good spacing

The key parameter of a real-space calculation is the spacing between the points of the mesh. The default option in {{< octopus >}} is to use a Cartesian grid. This is a regular grid where the spacing along each Cartesian direction is the same. The first step in any calculation should then be making sure that this spacing is good enough for our purposes. This should be done through a convergence study, very similar to the ones performed in plane-wave calculations.

The needed spacing essentially depends on the pseudopotentials that are being used. The idea is to repeat a series of ground-state calculations, with identical input files except for the grid spacing. There are many different ways of doing it, the simplest one being to change the input file by hand and run {{< octopus >}} each time. Another option is to use a little bash script like the one bellow. 

As a first example, we will use the Nitrogen atom from the previous tutorial, so make sure you have the corresponding {{< file "inp" >}} and {{< file "N.xyz" >}} files in your working directory. Next, create a file called {{< code spacing.sh >}} with the following content:


```bash
#include_file doc/tutorials/octopus_basics/total_energy_convergence/1.nitrogen_atom_spacing/spacing.sh
```

What this script does is quite simple: it runs {{< octopus >}} for a list of spacings (0.26, 0.24, ..., 0.14) and after each calculation it collects the values of the total energy and of the eigenvalues, and prints them to a file called {{< file "spacing.log" >}}. In order to change the value of the spacing between each calculation, it uses the feature that one can override variables in the {{< file "inp" >}} file by {{< manual "Advanced_ways_of_running_Octopus#passing-arguments-from-environment-variables" "defining them as environment variables" >}} in the shell. Note that we unset the variable OCT_Spacing at the end of the script, to avoid it to remaining defined in case you run the script directly in the shell. 

Now, to run the script type 
```bash
source spacing.sh
```
(**NOT** {{< command "sh spacing.sh" >}} or {{< command "bash spacing.sh" >}}!). 

Note that for this to work, the {{< file "octopus" >}} executable must be found in the shell path. If that is not the case, you can replace
```bash
  octopus >& out-$Spacing
```
with the actual path to the executable
```bash
  /path/to/octopus >& out-$Spacing
```

Once the script finishes running, the {{< file "spacing.log" >}} file should look something like this: 

```text
#include_file doc/tutorials/octopus_basics/total_energy_convergence/1.nitrogen_atom_spacing/spacing.log
```

You can also plot the results from the file in your favorite plotting program, e.g. {{< command "gnuplot" >}}.

#include_eps doc/tutorials/octopus_basics/total_energy_convergence/1.nitrogen_atom_spacing/NitrogenSpacing.eps caption="Convergence with spacing of N"

{{%expand "Expand to see the gnuplot script"%}}
```bash
#include_file doc/tutorials/octopus_basics/total_energy_convergence/1.nitrogen_atom_spacing/plot.gp
```
{{%/expand%}}

The results, for this particular example, are shown in the figure. In this figure we are actually plotting the error with respect to the most accurate results (smallest spacing). That means that we are plotting the '''difference''' between the values for a given spacing and the values for a spacing of 0.14 Å. So, in reading it, note that the most accurate results are at the left (smallest spacing). A rather good spacing for this nitrogen pseudopotential seems to be 0.18 Å. However, as we are usually not interested in total energies, but in energy differences, probably a larger one may also be used without compromising the results.

## Methane molecule

We will now move on to a slightly more complex system, the methane molecule CH<sub>4</sub>, and add a convergence study with respect to the box size.

### Input

As usual, the first thing to do is create an input file for this system. From our basic chemistry class we know that methane has a tetrahedral structure. The only other thing required to define the geometry is the bond length between the carbon and the hydrogen atoms. If we put the carbon atom at the origin, the hydrogen atoms have the coordinates given in the following input file:

{{< code-block>}}
#include_input doc/tutorials/octopus_basics/total_energy_convergence/2.methane/inp
{{< /code-block >}}

Here we define a variable CH that represents the bond length between the carbon and the hydrogen atoms, which simplifies the writing of the coordinates. We start with a bond length of 1.2 Å but that's not so important at the moment, since we can optimize it later (for more information see the {{<tutorial "Geometry optimization">}} tutorial). (You should not use 5 Å or so, but something slightly bigger than 1 Å is fine.)

Some notes concerning the input file:
* We do not tell {{< octopus >}} explicitly which {{< variable "BoxShape" >}} to use. The default is a union of spheres centered around each atom. This turns out to be the most economical choice in almost all cases.
* We also do not specify the %{{< variable "Species" >}} block. In this way, {{< octopus >}} will use default pseudopotentials for both Carbon and Hydrogen. This should be OK in many cases, but, as a rule of thumb, you should do careful testing before using any pseudopotential for serious calculations.

If you use the given input file you should find the following values in the resulting {{< file "static/info" >}} file.

{{< code-block >}}
#include_file doc/tutorials/octopus_basics/total_energy_convergence/2.methane/info.txt
{{< /code-block >}}


Now the question is whether these values are converged or not. This will depend on two things: the {{< variable "Spacing" >}}, as seen above, and the {{< variable "Radius" >}}. Just like for the Nitrogen atom, the only way to answer this question is to try other values for these variables.

###  Convergence with the spacing  


As before, we will keep all entries in the input file fixed except for the spacing that we will make smaller by 0.02 Å all the way down to 0.1 Å. So you have to run {{< octopus >}} several times. You can use a similar script as the one from the Nitrogen atom example:

```bash
#include_file doc/tutorials/octopus_basics/total_energy_convergence/3.methane_spacing/spacing.sh
```

For this example we changed the list of spacings to try and, to keep things simple, we are only writing the total energy to the {{< file "spacing.log" >}} file.

Now run the script to get the following results:

{{< code-block >}}
#include_file doc/tutorials/octopus_basics/total_energy_convergence/3.methane_spacing/spacing.log
{{< /code-block >}}


If you give these numbers to ''gnuplot'' (or other software to plot) you will get a curve like the one shown below:

#include_eps doc/tutorials/octopus_basics/total_energy_convergence/3.methane_spacing/Spacing_CH4.eps caption="Convergence with spacing of methane"


As you can see from this picture, the total energy is converged to within 0.1 eV for a spacing of 0.18 Å. So we will use this spacing for the next calculations.

### Convergence with the radius


Now we will see how the total energy changes with the {{< variable "Radius" >}} of the box. We will change the radius in steps of 0.5 Å. You can change the input file by hand and run Octopus each time, or again use a small script:

```bash
#include_file doc/tutorials/octopus_basics/total_energy_convergence/4.methane_radius/radius.sh
```


Before running the script, make sure that the spacing is set to 0.18 Å in the input file. You should then get a file {{< file "radius.log" >}} that looks like this:

```text
#include_file doc/tutorials/octopus_basics/total_energy_convergence/4.methane_radius/radius.log
```

Below you can see how this looks like when plotted. Note that in this case, the most accurate results are on the right (larger radius).

#include_eps doc/tutorials/octopus_basics/total_energy_convergence/4.methane_radius/Radius_CH4.eps caption="Convergence with radius of methane"

If we again ask for a convergence up to 0.1 eV we should use a radius of 3.5 Å. How does this compare to the size of the molecule? Can you explain why the energy increases (becomes less negative) when one decreases the size of the box?

{{< tutorial-footer >}}

