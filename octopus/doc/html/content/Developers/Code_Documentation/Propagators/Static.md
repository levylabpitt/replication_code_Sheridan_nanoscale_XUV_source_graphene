---
Title: "Static Propagator"
Weight: 9
---

The "static propagator" is just a dummy propagator, which does not propagate the system. The only implemented operation is to 
update the interactions, as this is necessary in a multi-system calculation, where other systems are propagated by a "real" propagator.


```Fortran
#include_type_def propagator_static_t

```
These are used to define the algorithm, which is done in the constructor of the propagator:
```Fortran
#include_function propagator_static_constructor
```


