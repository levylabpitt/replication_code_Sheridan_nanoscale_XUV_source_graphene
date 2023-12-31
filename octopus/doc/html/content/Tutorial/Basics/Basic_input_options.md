---
title: "Basic input options"
series: "Tutorials"
#tags: ["Basic", "Ground State", "Molecule", "Pseudopotentials", "DFT", "Total Energy"]
weight: 2
tutorials: ["Octopus Basics"]
difficulties: "basic"
difficulties_weight: 1
theories: "DFT"
calculation_modes: "Ground state"
system_types: "Molecule"
species_types: "Pseudopotentials"
features: "Total energy"
description: "Obtain the ground state of the nitrogen atom."
---


Now we will move to a more complicated (and realistic) input file. We will obtain the ground state of the nitrogen atom. We will introduce several basic input variables and will give a more detailed description of the output for this example.

## The input files

This sample input file lets us obtain the ground state of the nitrogen atom, within the LDA approximation, in a closed-shell (unpolarized) configuration (as explained below, you need an auxiliary {{< file ".xyz" >}} input). Note that this is not the correct ground state of the nitrogen atom! However, it will permit us to describe some of the most important input variables:

```
#include_input testsuite/tutorials/02-octopus_basics-basic_input_options.01-N_atom.inp
```

We have introduced here several new variables:

* {{< code-inline >}}{{< variable "UnitsOutput" >}} = eV_Angstrom{{< /code-inline >}}: Two different unit systems may be used for output: the usual atomic units (which is the default, and the ones used internally in the code); and the system in which the Ångström is substituted for the atomic unit of length, and the electronvolt is substituted for the atomic unit of energy. You can find a more detailed description of units in {{< octopus >}} in the {{< manual "Basics/Units" "Units" >}} page of the manual.

* The following entry in the input file is not a variable that {{< octopus >}} will read directly, but rather illustrates the possibility of writing "user-defined" values and expressions to simplify the input file. In this case, we define the nitrogen mass ({{< code-inline >}}Nitrogen_mass = 14.0{{< /code-inline >}}) (note that in this case, as an exception, the value is expected to be in the so-called "atomic mass units", rather than in "atomic units"). This definition may be used elsewhere in the input file.

* The {{< variable "Species" >}} block should contain the list of species that are present in the system to be studied. In this case we have only one species: nitrogen. The first field is a string that defines the name of the species, "N" in this case. The second field defines the type of species, in this case {{< code-inline >}}species_pseudo{{< /code-inline >}}. Then a list of parameters follows. The parameters are specified by a first field with the parameter name and the field that follows with the value of the parameter. Some parameters are specific to a certain species while others are accepted by all species. In our example {{< code-inline >}}set{{< /code-inline >}} instructs {{< octopus >}} to use a pseudopotential for nitrogen from the {{< code-inline >}}standard{{< /code-inline >}} set. This happens to be a Troullier-Martins pseudopotential defined in the {{< file "N.psf" >}} file found in the directory {{< file "share/octopus/pseudopotentials/PSF" >}}. Then come maximum {{< code-inline >}}lmax{{< /code-inline >}} - component of the pseudopotential to consider in the calculation, and the {{< code-inline >}}lloc{{< /code-inline >}} - component to consider as local. Generally, you want to set the maximum ''l'' to the highest available in the pseudopotential and the local ''l'' equal to the maximum ''l''. Finally, the mass of the species can also be modified from the default values by setting {{< code-inline >}}mass{{< /code-inline >}} parameter.

* {{< code-inline >}}{{< variable "XYZCoordinates" >}} = 'N.xyz'{{< /code-inline >}}: The geometry of the molecule (in this case, a single atom in the grid origin) is described in this case in a file with the well known {{< code-inline >}}XYZ{{< /code-inline >}} format. The file for this outrageously simple case is given by:

{{< code-block >}}
#include_input doc/tutorials/octopus_basics/basic_input_options/N.xyz
{{< /code-block >}}

* {{< code-inline >}}{{< variable "ExtraStates" >}} = 1{{< /code-inline >}}: By default, {{< octopus >}} performs spin-unpolarized calculations (restricted closed-shell, in Hartree-Fock terminology). It then places two electrons in each orbital. The number of orbitals, or Kohn-Sham states, is then calculated by counting the number of valence electrons present in the system, and dividing by two. In this case, since we have five valence electrons, the code would use three orbitals. However, we know beforehand that the HOMO orbital has a three-fold degeneracy, and as a consequence we need to put each one of the three _p_ electrons in a different orbital. We therefore need one more orbital, which we get with this line in the input file.

* {{< code-inline >}}%{{< variable "Occupations" >}}{{< /code-inline >}} block: Generally, the occupations of the Kohn-Sham orbitals are automatically decided by the code, filling the lowest-energy orbitals. However, if we have degeneracies in the LUMO as in this case, the user may want to accommodate the electrons in a certain predefined way. In this example, the obvious way to fill the orbitals of the nitrogen atom is to put two electrons in the first and deepest orbital (the _s_ orbital), and then one electron on each of the second, third and fourth orbitals (the _p_ orbitals, which should be degenerate).

* {{< code-inline >}}{{< variable "BoxShape" >}} = sphere{{< /code-inline >}}: This is the choice of the shape of the simulation box, which in this case is set to be a sphere (other possible choices are {{< code-inline >}}minimum{{< /code-inline >}}, {{< code-inline >}}cylinder{{< /code-inline >}}, or {{< code-inline >}}parallelepiped{{< /code-inline >}}).


## Output

Once you have constructed the input file and created the {{< file "N.xyz" >}} file, you may unleash {{< octopus >}} on it. Lets now go over some of the sections of the output.

#### Species

{{< code-block >}}
#include_file doc/tutorials/octopus_basics/basic_input_options/Species.txt
{{< /code-block >}}

Here the code searches for the needed pseudopotential files, and informs the user about its success or failure. In this case, only the {{< file "N.psf" >}} file is required. Once that file has been processed, some information about it is written to the output. One of the most important pieces of information to be found here is the valence charge, which tells us how many electrons from this species will be considered in the calculation.

#### Grid

{{< code-block >}}
#include_file doc/tutorials/octopus_basics/basic_input_options/Grid.txt
{{< /code-block >}}

This step is about the construction of the mesh. As requested in the input file, a sphere of radius 5 Å is used, which contains a cubic regular real-space grid with spacing 0.18 Å. This implies 89727 points ({{< code-inline >}}inner mesh =  89727{{< /code-inline >}}). For the sake of comparison with plane-wave-based codes, this is more or less equivalent to a plane-wave calculation that imposes a density cutoff of 1160.595 eV = 42.6 Hartree (except that in this case there is no artificial periodic repetition of the system).

#### Mixing

{{< code-block >}}
#include_file doc/tutorials/octopus_basics/basic_input_options/Mixing.txt
{{< /code-block >}}

During the self-consistent procedure one has to use a {{< manual "Calculations/Ground_State#Mixing" "mixing scheme" >}} to help convergence. One can mix either the density or the potential, and there are several mixing schemes available.

#### Eigensolver

{{< code-block >}}
#include_file doc/tutorials/octopus_basics/basic_input_options/Eigensolver.txt
{{< /code-block >}}

Here we see that the {{< manual "Calculations:Ground_State#Eigensolver" "eigensolver" >}} used will be simple conjugate gradients (cg), and a preconditioner is used to speed up its convergence.

#### LCAO
After some output you should see something like:

{{< code-block >}}
#include_file doc/tutorials/octopus_basics/basic_input_options/lcao.txt
{{< /code-block >}}

This is the first step of a ground-state calculation: obtaining a reasonably good starting density and Kohn-Sham orbitals to feed in the self-consistent (SCF) procedure. For this purpose, {{< octopus >}} performs an initial calculation restricted to the basis set of atomic orbitals ( {{< manual "Calculations:Ground_State#LCAO" "Linear Combination of Atomic Orbitals" >}}, LCAO). The resulting eigenvalues of this calculation are written to standard output.

#### Wavefunction kind

{{< code-block >}}
 Info: SCF using real wavefunctions.
{{< /code-block >}}

Very often one can work with real wave-functions. This is particularly helpful as calculations with real wave-functions are much faster than with complex ones. However, if a magnetic field is present, if the system is periodic, or if spin-orbit coupling is present, complex wave-functions are mandatory. But don't worry: the program is able to figure out by itself what to use.

#### SCF

{{< code-block >}}
#include_file doc/tutorials/octopus_basics/basic_input_options/first_iter.txt
{{< /code-block >}}


Now the SCF cycle starts. For every step, {{< octopus >}} outputs several pieces of information:

* The values {{< code-inline >}}abs_dens{{< /code-inline >}} and {{< code-inline >}}rel_dens{{< /code-inline >}} are to monitor the absolute and relative convergence of the density, while {{< code-inline >}}rel_ev{{< /code-inline >}} and {{< code-inline >}}abs_ev{{< /code-inline >}} are two alternative measures of the convergence, based on measuring the difference between input and output eigenvalues. The SCF procedure, by default, is stopped when {{< code-inline >}}rel_dens{{< /code-inline >}} is smaller than $10^{-5}$. This may be altered with the appropriate input variables (see in the manual the variables {{< variable "ConvAbsDens" >}}, {{< variable "ConvRelDens" >}}, {{< variable "ConvAbsEv" >}} and {{< variable "ConvRelEv" >}}).

* The line {{< code-inline >}}Matrix vector products:    108{{< /code-inline >}} tells us that the Hamiltonian was applied 108 times. This gives us an idea of the computational cost.

* The line {{< code-inline >}}Converged eigenvectors:      0{{< /code-inline >}} tells us that upon completion of the diagonalization procedure, none of the orbitals met the required precision criterion for the wavefunctions. In a following example, we will modify this criterion in the input file.

* The list of eigenvalues is then printed, along with their errors: how much they deviate from "exact" eigenvalues of the current Hamiltonian. This number is the so-called "residue".

You can now take a look at the file {{< file "static/info" >}} that will hold a summary of the calculation.

## Restarting

Any ground-state calculation may be restarted later (to refine it if it did not converge properly, or with any other purpose), provided that the contents of the {{< code restart >}} directory are preserved. You can try this now, just by running {{< octopus >}} again. You will notice that {{< octopus >}} did not give any warning after the line

{{< code-block >}}
 Info: Loading restart information.
{{< /code-block >}}

This is useful if you change slightly the parameters of the simulation (for example the XC functional or the convergence criteria). If you change the grid parameters {{< octopus >}} will not be able to restart from the previous calculation. If you do not want {{< octopus >}} to try to restart a calculation, you can set the variable {{< variable "FromScratch" >}}.

In case you ware wondering what the restart information looks like, you can have a look at the contents of the {{< file "restart" >}} directory. This is where the files needed to restart a calculation are stored. It may contain several sub-directories depending on the calculations previously performed. In this case, it just contains one:

```bash
 % ls restart
 gs
 % ls restart/gs
#include_file doc/tutorials/octopus_basics/basic_input_options/ls_restart_gs.txt
```

{{< octopus >}} stores each individual state in a different binary (yet platform-independent) file. In this case, we only have four states (files {{< file "0000000001.obf" >}} to {{< file "0000000004.obf" >}}). Some other useful quantities, like the density, are also stored in binary form. The other files are text files that contain diverse control information. It is unlikely that you will ever have to work directly with these files, but you may take a look around if you are curious. 


{{< tutorial-footer >}}

