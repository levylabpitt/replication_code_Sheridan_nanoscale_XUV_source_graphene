---
Title: "Calculating energies"
Weight: 15
---

In the multisysten framework, the total energy ({{<code "total_energy">}}) of a calculation consists of various contributions, where we follow the standard definitions of
thermodynamics (see [here](https://en.wikipedia.org/wiki/Internal_energy))

* The kinetic energy ({{<code "kinetic_energy">}})
* The internal energy ({{<code "internal_energy">}})
* The potential energy ({{<code "potential_energy">}})

### Kinetic energy

We use the term *kinetic energy* in the usual sense. The only slight exception is for container systems:
{{% notice note %}}
Ideally, for containers, the kinetic energy should be the kinetic energy with respect to the centre of mass of the container.
The kinetic energies of the constituents, relative to the centre of mass frame, should be accounted for as part of the internal energy of the container.
{{% /notice %}}

{{% expand "multisystem_update_kinetic_energy()" %}}
```Fortran
#include_subroutine multisystem_update_kinetic_energy
```
{{% /expand %}}

Specific systems need their specific routines to calculate the kinetic energy.
{{% notice note %}}
For the Maxwell system, this is the so-called electromagnetic energy.
{{% /notice %}}

### Interaction energies

Everything else is treated as an interaction energy. Here we have to distinguish between the interactions between two different systems (which can be of the same type), and the interaction of a system with itself (intra-interaction). 

#### inter-interactions

The interactions between two different systems do not require any further thought, as by definition no physical self-interaction can occur. 
As the method to calculate the interactions are part of the interaction class, it is independent of the actual system and can be implemented 
at the level of the class{{<code "system_t">}}.

{{% expand "system_update_potential_energy()" %}}
```Fortran
#include_subroutine system_update_potential_energy
```
{{% /expand %}}

The exception are containers. Here we need to loop over the constituents. In order to distinguish inter- from intra-interactions, we need to
query each interaction for its interaction partner, and skip the interaction, if the partner is part of the container.

{{% expand "multisystem_update_potential_energy()" %}}
```Fortran
#include_subroutine multisystem_update_potential_energy
```
{{% /expand %}}

#### intra-interactions

Systems may contain more than one physical particle (e.g. the electrons, a set of ions or container systems). In order to account for the interaction of these particles with other particles of the same system, we decided to treat this case as a system interacting with itself, which we call intra-interaction. 

In some cases, such as a single particle, this intra interaction has to be zero, while in other cases with many particles, the interactions have to be calculated, where -- of course -- the interaction of one particle with itself has to be removed (at least approximatively).


Another important aspect of the implementation is that {{<octopus>}} deals with _one-sided_ interactions. This has the implication that there is no double counting when calculating
the interaction energies. Both contributions have to be counted: 
**for each system, we add the the energy of the system in the field of the partner.**

{{% expand "system_update_internal_energy()" %}}
```Fortran
#include_subroutine system_update_internal_energy
```
{{% /expand %}}


{{% notice note %}}
For containers, this means that the interaction energy contains the complete interaction energies of all subsystems in the container, _plus_ the energies of the subsystems in the field of all other systems (not part of that container).
{{% /notice %}}

{{% notice warning %}}
Note that the internal_energy of a container is _not_ the sum of the self_energies of the constituents, and also the external energy of a container is _not_
the sum of the external energies of the constituents!
{{% /notice %}}

{{% expand "multisystem_update_internal_energy()" %}}
```Fortran
#include_subroutine multisystem_update_internal_energy
```
{{% /expand %}}





### Total energy

The total energy of the whole simulation (top level container) is clearly defined. It contains the sum of all internal and interaction energies.

For a specific system (e.g. electrons), the total energy also is the sum of the internal energy and the interaction energies (including the intra interaction energy).

{{% expand "system_update_total_energy()" %}}
```Fortran
#include_subroutine system_update_total_energy
```
{{% /expand %}}

 

{{% notice warning %}}
Total energies and interaction energies depend on the state of other systems, and therefore should not be used for convergence criteria for a given system. Only the internal energy is safe to use here.
{{% /notice %}}

