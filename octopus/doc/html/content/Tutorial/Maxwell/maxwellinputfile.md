---
Title: "Maxwell input file"
series: "Tutorials"
tutorials: "Maxwell"
author: ["Franco Bonafé","René Jestaedt","Heiko Appel","Martin Lueders"]
weight: 4
---

### Maxwell Input File

#### Input file variable description

##### Calculation mode and parallelization strategy


At the beginning of the input file, the basic variable {{< variable
"CalculationMode" >}} selects the run mode of Octopus and has to be set always.
In case of a parallel run, there are some variables to set the proper
parallelization options. For Maxwell propagations, parallelization in domains
is possible, but parallelization in states is not needed, as there are always 3
or 6 states, depending on the Hamiltonian (see below).

{{< code-block >}}
#include_input_snippet doc/tutorials/maxwell/00-maxwell-input-options.inp calc_mode
{{< /code-block >}}


##### Multisystem setup

The Maxwell system is implemented as part of the Octopus multisystem framework.
This allows to calculate several different systems, which can interact with
each other.

Currently implemented system types are:

* electronic: An electronic system. (only partly implemented)
* maxwell: A maxwell system.
* classical_particle: A classical particle. Used for testing purposes only.
* charged_particle: A charged classical particle.
* multisystem: A container system containing other systems.
* dftbplus: A DFTB+ system
* linear_medium: A (static) linear medium for classical electrodynamics.
* dispersive_medium: A (dispersive) linear medium (only the Drude model is currently implemented)


## Definition of the systems

In this tutorial, we will use a pure Maxwell system:

{{< code-block >}}
#include_input_snippet doc/tutorials/maxwell/00-maxwell-input-options.inp system
{{< /code-block >}}

Subsequently, variables relating to a named system are prefixed by the system
name. However, this prefix is not mandatory for calculations considering only
one system.

## Maxwell box variables and parameters

### Simulation box

The Maxwell simulation box is of parallelepiped shape, and is divided into mainly two
regions: one inner region for the requested propagation and an outer region to simulate the
proper boundary conditions. In the input file, the simulation box sizes
refer to the total simulation box which includes the boundary region. The
inner simulation box is defined by the total box size minus the boundary
width. In case of zero boundary condition, there is no boundary region. {{<
figure caption="Boundaries in 2D" src="/images/Maxwell/boundaries_in_2D.png"
width="500px" >}}

The boundary region can be set up by absorbing boundaries, incident plane waves
or a combination of both. In the latter case, the absorbing boundary region is
the inner boundary region, and the plane waves region is the outer boundary
region. Again, the box sizes are determined by the total simulation box size
and the corresponding boundary width. {{< figure caption="Plane Waves and PML"
src="/images/Maxwell/plane_wave_and_pml_boundaries_in_2D.png" width="500px" >}}

{{< notice note >}} For a given size of the propagation region, the boundaries
have to be added to the total box size. The boundary width is given by the
derivative order (default is 4) times the spacing. The width of the absorbing
boundary region is defined by the variable {{< code "MaxwellABWidth" >}}.
{{< /notice >}}

<!---
The matter grid is in general located inside the Maxwell grid. There are
several possible types of grids to describe a coupled Maxwell-matter system.
The following figure illustrates some possible cases overlaying Maxwell and matter
grids. In the Octopus code, currently only the types e), f), and g) are
implemented. So the matter box sizes and Maxwell box sizes can be chosen
independently, whereas the spacings of both grids have to be equal and the grid
points have to lie on the top of each other. The only exception is type g),
where the matter grid is much smaller than the Maxwell grid. In this case, the
matter grid size has to be smaller than the Maxwell grid spacing. {{< figure
caption="Multiscale" src="/images/Maxwell/multiscale_figures.png" width="500px"
>}}
--->


##### Maxwell options

First, the Maxwell propagation operator options consist of the type of {{< code
"MaxwellHamiltonianOperator" >}}, which normally will be a
Maxwell Hamiltonian for propagation in vacuum (faraday_ampere). This options is also used
for dispersive media. For linear (static) media, one should set the Maxwell-Hamiltonian
that works for a linear medium (faraday_ampere_medium).

Currently, the Maxwell system does not define its own propagator, but uses the
system default propagator, defined in {{< code "TDSystemPropagator" >}}. So
far, only the exponential midpoint is implemented.

The Maxwell Boundary conditions can be chosen for each direction differently
with {{< code "MaxwellBoundaryConditions" >}}. Possible choices are zero,
plane waves and constant field boundary conditions. Additionally to the boundary condition,
the code can include absorbing boundaries. This means that all outgoing waves are absorbed
while all incoming signals still arise at the boundaries and propagate into
the simulation box. The absorbing boundaries can be achieved by a mask
function or by a perfectly matched layer (PML) calculation with additional
parameters. For more information on the physical meaning of these parameters,
please refer to the [Maxwell-TDDFT paper].

[Maxwell-TDDFT paper]: https://doi.org/10.1080/00018732.2019.1695875

{{< code-block >}}
#include_input_snippet doc/tutorials/maxwell/00-maxwell-input-options.inp maxwell_calc
{{< /code-block >}}

##### Output options

The output option {{< code "OutputFormat" >}} can be used exactly as in the
case of TDDFT calculations, with the exception that some formats are inconsistent
with a Maxwell calculation (like xyz). It is possible to chose multiple
formats. The Maxwell output options are defined in the {{< code
"MaxwellOutput" >}} block (written to output_iter every {{< code
"MaxwellOutputInterval" >}} steps) and through the {{< code
"MaxwellTDOutput" >}} (written to td.general every step).

{{< code-block >}}
#include_input_snippet doc/tutorials/maxwell/00-maxwell-input-options.inp output2
{{< /code-block >}}


##### Time step variables

The {{< code "TDTimeStep" >}} option in general defines the time step used
for each system. For the stability of the time propagation,
Maxwell systems need to fulfill the Courant condition. In some cases larger
time steps are possible but this should be checked case by case.

{{< code-block >}}
#include_input_snippet doc/tutorials/maxwell/00-maxwell-input-options.inp timestep
{{< /code-block >}}


##### Maxwell field variables

The incident plane waves can be defined at the boundaries using the {{<
code "MaxwellIncidentWaves" >}} block. This block defines the incoming wave
type and complex electric field amplitude, and an envelope function
name, which must be defined in a separate {{<
code "MaxwellFunctions" >}} block. Multiple plane waves can be
defined, and their superposition gives the final incident wave result at the
boundary. This plane waves boundaries are evaluated at each time step, and give
the boundary conditions for the propagation inside the box. To use them {{<
code "MaxwellBoundaryConditions" >}} must be set to "plane_waves" in the
proper directions.

The {{< code "MaxwellFunctions" >}} block reads the additional required
information for the plane wave pulse, for each envelope function name given in
the {{< code "MaxwellIncidentWaves" >}} block.

Inside the simulation box, the initial electromagnetic field is set up by the
{{< code "UserDefinedInitialMaxwellStates" >}} block. If the pulse
parameters are such that it is completely outside the box at the initial time,
this variables does not need to be set, as the default initial field inside the
box is zero.

{{< notice warning >}}
If the incident wave pulse should have initially a non-zero value
inside the box, and {{< code "UserDefinedInitialMaxwellStates" >}} is not set,
the initial field will be zero by default. This means having a discontinuity
at the boundaries for t=0, which will lead to an erroneous result!
{{< /notice >}}

In case of {{< code-inline >}} MaxwellPlaneWavesInBox = yes{{< /code-inline >}},
no Maxwell propagation is done, but instead the field
is set analytically at each time step using the parameters defined in {{<
code "MaxwellIncidentWaves" >}}.

{{< code-block >}}
#include_input_snippet doc/tutorials/maxwell/00-maxwell-input-options.inp field
{{< /code-block >}}

##### External current

An external current can be added to the internal current. Its spatial and
temporal profile are analytical, and the formulas are defined by the
{{< code "UserDefinedMaxwellExternalCurrent" >}} block. This block defines
the spatial profile and frequency of each current source, and it sets an
envelope function name, that must be defined in a separate {{< code
"TDFunctions" >}} block.

{{< code-block >}}
#include_input_snippet doc/tutorials/maxwell/00-maxwell-input-options.inp ext_current
{{< /code-block >}}


##### Linear Medium

Different shapes of dispersive and non-dispersive linear media can be included
in the calculation either through a pre-defined box shape (parallelepiped,
sphere or cylinder), or through a file describing more complex geometries. For
more information check the Linear Medium tutorial.


{{< tutorial-footer >}}

