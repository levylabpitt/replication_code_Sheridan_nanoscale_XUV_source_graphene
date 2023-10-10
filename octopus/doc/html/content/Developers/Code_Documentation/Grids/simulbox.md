---
Title: Simulation Box
Weight: 1
---


The simulation box is defined by the class:
```Fortran
#include_type_def simul_box_t
```
which is derived from {{< developers "Code_documentation/ions/box" "box_t" >}}.

The main function of the simulation box is {{< code "simul_box_contains_points" >}}, which determines for a set of points whether they are inside or outside the simulation box.

{{% expand "Definition of simul_box_contains_points()" %}}
```Fortran
#include_function simul_box_contains_points
```
{{% /expand %}}
