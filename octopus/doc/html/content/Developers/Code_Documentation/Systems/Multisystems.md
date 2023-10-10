---
Title: "Multisystem: Containers"
Weight: 2
---

### Introduction

In {{<octopus>}} there is the possibility to group several systems into containers, which are implemented in the {{<code "multisystem_basic_t">}} class.

Containers can have different purposes and applications. In the simplest case, containers are simply a collection of other systems, and do not have their own interactions with anything else. In this case, containers do not introduce any different physics (or approximations), but simply help in the book-keeping of the problem.

{{% notice note %}}
Containers do not correspond to a given region in space, but only to a selection of systems. In many cases, these systems *might* be confined to a certain region in space, but this is not a property of the container. In many other cases, e.g. combining matter and maxwell fields, both systems occupy the same space, or have a substantial overlap.
{{% /notice %}}

Another use case might be to group systems into a container, and then only interact with the whole container, instead of the individual systems. This, however, is an approximation, and furthermore (at least, at the moment) has some limitations due to the implementation.

{{% notice warning %}}
Note, that containers themselves do not move. This has consequences to the definition of the energy contributions. In particular, a container does not have it's own kinetic energy, and all kinetic energy contributions of the constituents are accounted for in the internal energy. For more information, see 
{{<versioned-link "Developers/Code_Documentation/Multisystem_framework/energies" "Calculating Energies">}}
{{% /notice %}}


### Implementation

Most of the functionality is implemented at the level of the abstract class {{<code "multisystem_t">}}:
{{% expand "Definition of multisystem_t"%}}
```Fortran
#include_type_def multisystem_t
```
{{% /expand %}}

The specific {{<code "multisystem_basic_t">}} class only adds the finalizer:
{{% expand "Definition of multisystem_basic_t"%}}
```Fortran
#include_type_def multisystem_basic_t
```
{{% /expand %}}
