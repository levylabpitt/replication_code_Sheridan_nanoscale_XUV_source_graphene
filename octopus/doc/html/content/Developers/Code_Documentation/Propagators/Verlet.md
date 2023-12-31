---
Title: "Detailled example: Verlet"
Weight: 10
---

{{< notice warning >}}
Work in progress!
{{< /notice >}}

### The Verlet algorithm

According to [Wikipedia](https://en.wikipedia.org/wiki/Verlet_integration), the Verlet algorithm is defined as:

- Calculate $\vec{x}(t + \Delta t) = \vec{x}(t) + \vec{v}(t) \Delta t + \tfrac12 \vec{a}(t) \Delta t^2$.
- Derive $\vec{a}(t + \Delta t)$ from the interaction potential using $\vec{x}(t + \Delta t)$.
- Calculate $\vec{v}(t + \Delta t) = \vec{v}(t) + \tfrac12 \big(\vec{a}(t) + \vec{a}(t + \Delta t)\big)\Delta t$.



The Verlet operator is represented by the type {{< code "propagator_verlet_t" >}} which extends {{< code "propagator_t" >}}:

{{% expand "Definition of propagator_verlet_t" %}}
```Fortran
#include_type_def propagator_verlet_t
```
{{% /expand %}}

{{< notice note >}}
Note, that the class definition does not add anything to the {{< code "propagator_t" >}} class. The only differences are the definition of the operations, and the overloaded constructor. 
{{< /notice >}}

For the Verlet propagator, we need to define the following operations:
```Fortran
#include_code_doc verlet_propagation_operations
```
These are used to define the algorithm, which is done in the constructor of the propagator:
```Fortran
#include_function propagator_verlet_constructor
```
The algorithm also uses steps, which are not specific to the Verlet algorithm, and are defined in the {{< code "propagator_oct_m" >}} module.

{{< notice note >}}
Note, the difference between {{< code "OP_VERLET_FINISH" >}} and {{< code "OP_FINISHED" >}}. The former denotes a specific step, to be taken at the end of one time step, 
while the latter generally denotes the end of a time step.
{{< /notice >}}

As can be seen in this example, the definition of the propagator in terms of the algorithmic operations is a direct translation of the algorithm, shown above.

### Implementation of the steps

So far, the propagator is only defined in an abstract way. It is down to the individual systems (or better, their programmers) to implement the individual steps of the algorithm, in terms of the dynamic variables for that system. An easy examample to demonstrate this is the classical particle, as implemented in {{< code "classical_particles_t" >}}

{{% expand "definition of classical_particles_t" %}}
```Fortran
#include_type_def classical_particles_t
```
{{% /expand %}}
This describes an extremely simple system, consisting of a set of classical, neutral particles. The dynamic variables (i.e. the state) of the system are the positions, velocities and accelerations.

One ''tick'' of the propagator is defined in the function {{< code "classical_particle_do_td" >}}. As can be seen in the code below, this function implements all possible algorithmic steps for all propagators, allowed for that system.


{{% expand "Example implementation for classical particles" %}}
```Fortran
#include_subroutine classical_particles_do_td
```
{{% /expand %}}


### The timeline explained

The Verlet algorithm is a good (because simple) example to illustrate how the systems and the interaction are updated as the code progresses through the main loop.


{{< d3-sequence file="graph_data/propagation-3body-verlet-equal-step.json"  viewContainers="yes" viewGhosts="yes" >}}

This graph illustrates how the state machine is stepping through the algorithm. Each system is picking the next algorithmic step from the propagator. For the containers (i.e. ''root'' and ''earth''), the only steps are ''Updating interactions'' and ''Finished''. The {{< emph real >}} systems, on the other hand, are progressing orderly through the operations, defined in the propagator.


</br>

{{% expand "Example with different time steps" %}}
This example shows the propagation of the solar system, where different time steps are used for the three systems.

{{< d3-sequence file="graph_data/propagation-3body-verlet-different-step.json" viewContainers="yes" viewGhosts="yes" >}}

In contrast to the previous example we can see here that the slowest system (i.e. ''sun'') is waiting in ''Updating interaction'' for many computational steps, as it needs to wait for the other systems to complete their time steps. These waiting steps are computationally very cheap, as no grid operations are performed.

{{% /expand %}}
