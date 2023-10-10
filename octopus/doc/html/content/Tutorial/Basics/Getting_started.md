---
title: "Getting started"
authors: ["Micael Oliveira", "Xavier Andrade", "David Strubbe"]
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
weight: 1
description: "Learn how to run the code"
---


The objective of this tutorial is to give a basic idea of how {{< octopus >}} works.

### Generating the input file

With a text editor, create a text file called {{< file "inp" >}} containing the following text:

{{< code-block >}}
#include_input doc/tutorials/octopus_basics/getting_started/1.H_atom/inp
{{< /code-block >}}

This is the simplest example of an {{< octopus >}} input file:

* {{< code-inline >}}{{< variable "CalculationMode" >}} = gs {{< /code-inline >}}: This variable defines the run mode -- please consult the manual for the full list of the possible run modes. In this case we set it to {{< code gs >}}, which instructs the code to start a ground-state calculation.

* {{< code-inline >}}%{{< variable "Coordinates" >}}{{< /code-inline >}}: The entry is not just the definition of a variable, but rather of a full set of them -- a "block" of variables. The beginning of a block is marked by the {{< code "%identifier" >}} line, and ended by a {{< code "%" >}} line. In this case the identifier is {{< code-inline >}}%{{< variable "Coordinates" >}}{{< /code-inline >}}, where we list the atoms or species in our calculation and its coordinates, one per line. In this case, we put a single hydrogen atom in the center of our simulation box. 

* {{< code-inline >}}{{< variable "Spacing" >}} = 0.25*angstrom{{< /code-inline >}}: As you should know, {{< octopus >}} works in a real-space regular cubic mesh. This variable defines the spacing between points, a key numerical parameter, in some ways equivalent to the energy cutoff in plane-wave calculations.

* {{< code-inline >}}{{< variable "Radius" >}} = 4.0*angstrom{{< /code-inline >}}: The radius of the sphere that defines the simulation box.

The reason this input file can be so simple is that {{< octopus >}} comes with default values for the simulation parameters, and a set of default pseudopotentials for several elements (for properly converged calculations you might need to adjust these parameters, though).

To get a general idea of the format of the {{< octopus >}} input file, go and read the page about the {{< manual "Basics/Input file" "Input file" >}} in the manual.

The documentation for each input variable can be found in the {{< versioned-link "Variables/" "variable reference" >}} online, and can also be accessed via the {{< manual "utilities/oct-help" "oct-help" >}} utility.

### Running Octopus

Once you have written your input file, run the {{< command "octopus" >}} command (using {{< command "mpirun" >}} and perhaps a job script if you are using the parallel version). If everything goes correctly, you should see several lines of output in the terminal (if you don't, there must be a problem with your installation). As this is probably the first time you run {{< octopus >}}, we will examine the most important parts of the output.

{{% notice note %}}
Be aware that the precise values you find in the output might differ from the ones in the tutorial text. This can be due to updates in the code, or also changes in the compilation and run configuration.
{{% /notice %}}

* First there is an octopus drawn in ASCII art, the copyright notice and some information about the octopus version you are using and the system where you are running:
{{< code-block >}}
#include_file doc/tutorials/octopus_basics/getting_started/1.H_atom/header.txt
{{< /code-block >}}
Note that it also gives you the revision number, the compiler, and the compiler flags used. You should always include this information when submitting a bug report!

* The type of calculation it was asked to perform:
{{< code-block >}}
#include_file doc/tutorials/octopus_basics/getting_started/1.H_atom/Calculation_mode.txt
{{< /code-block >}}

* The spatial dimensions and the periodicity of the system:
{{< code-block >}}
#include_file doc/tutorials/octopus_basics/getting_started/1.H_atom/Space.txt
{{< /code-block >}}

* The species and pseudopotentials it is using:
{{< code-block >}}
#include_file doc/tutorials/octopus_basics/getting_started/1.H_atom/Species.txt
{{< /code-block >}}


* After some other output, {{< octopus >}} prints information about the grid: as we didn't say anything in the input file, {{< octopus >}} used the parameters recommended for this pseupopotential:
{{< code-block >}}
#include_file doc/tutorials/octopus_basics/getting_started/1.H_atom/Grid.txt
{{< /code-block >}}


* The level of theory and, in the case of (TD)DFT, the approximation to the exchange-correlation term:
{{< code-block >}}
#include_file doc/tutorials/octopus_basics/getting_started/1.H_atom/Theory_level.txt
{{< /code-block >}}


* At this point, {{< octopus >}} tries to read the wave-functions from a previous calculation. As there are none, it will give a warning.
{{< code-block >}}
** Warning:
**   Could not find 'restart/gs' directory for restart.
**   No restart information will be read.

** Warning:
**   Unable to read wavefunctions.
**   Starting from scratch!
{{< /code-block >}}

* Now {{< octopus >}} commences the calculation. To get a reasonable starting point for the DFT calculation, the initial wavefunctions are calculated as a {{< manual "Calculations:Ground_State#LCAO" "Linear Combination of Atomic Orbitals" >}} (LCAO).
{{< code-block >}}
#include_file doc/tutorials/octopus_basics/getting_started/1.H_atom/lcao.txt
{{< /code-block >}}

* After the LCAO, the real DFT calculation starts. For each self-consistency step some information is printed. When SCF {{< manual "Calculations:Ground_State#Convergence" "converges" >}}, the calculation is done.
{{< code-block >}}
#include_file doc/tutorials/octopus_basics/getting_started/1.H_atom/first_iter.txt
{{< /code-block >}}
...
{{< code-block >}}
#include_file doc/tutorials/octopus_basics/getting_started/1.H_atom/last_iter.txt
#include_file doc/tutorials/octopus_basics/getting_started/1.H_atom/footer.txt
{{< /code-block >}}

Just running the command {{< command "octopus" >}} will write the output directly to the terminal. To have a saved copy of the output, it is generally advisable to redirect the output into a file, and to capture the standard error stream as well, which can be done like this: {{< command "octopus &> log" >}}. That would create a file called {{< file "log" >}} containing all output including warnings and errors in their context.

###  Analyzing the results  

After finishing the calculation you will find a series of files in the directory you ran:

{{< code-block >}}
% ls
  exec inp restart static
{{< /code-block >}}

For the moment we will ignore the '''exec'''  and  '''restart''' directories and focus on the {{< file "static/info" >}} file, which contains the detailed results of the ground-state calculation. If you open that file, first you will see some parameters of the calculations (that we already got from the output) and then the calculated energies and eigenvalues in Hartrees:

{{< code-block >}}
#include_file doc/tutorials/octopus_basics/getting_started/1.H_atom/info.txt
{{< /code-block >}}


Since by default {{< octopus >}} does a spin-unpolarized density-functional-theory calculation with the local-density approximation, our results differ from the exact total energy of 0.5 H. Our exchange-correlation functional can be set by the variable {{< variable "XCFunctional" >}}, using the set provided by the {{% libxc %}} library.

### Extra

If you want to improve the LDA results, you can try to repeat the calculation with spin-polarization:

{{< code-block >}}
 {{< variable "SpinComponents" >}} = spin_polarized
{{< /code-block >}}

And if you want to obtain the exact Schödinger equation result (something possible only for very simple systems like this one) you have to remove the self-interaction error (a problem of the LDA). Since we only have one electron the simplest way to do it for this case is to use independent electrons:

{{< code-block >}}
 {{< variable "TheoryLevel" >}} = independent_particles
{{< /code-block >}}

A more general way would be to include self-interaction correction.

{{< tutorial-footer >}}

