---
Title: Mesh operations
Weight: 5
---



Accessing functions
===================

Function, e.g. the density are given by an array rho(_i_,_spin_), where _i_ ranges over the mesh points associated with the given domain (see domain decomposition), representing

rho( **r**<sub>_i_</sub>,  _spin_ )


What is `mesh_x_global(mesh, ip)` ?

`mesh_x_global(mesh, ip)` returns a `FLOAT` vector, corresponding to the coordinates of point `ip` of the global mesh.

```Fortran
#include_function mesh_x_global
```


Integration
-----------

Integration is implemented as a straightforward summation:

```Fortran
 do ip = 1, mesh%np
      dd = dd + ff(ip)
 end do 
```

Differentiation
---------------

Taking derivatives is done by finite differences. The points involved
in a derivative are defined by the stencil (see below).

Derivatives are discussed in a separate document [Derivatives.md](../Derivatives).


Note on packed states:
----------------------



Note on curvilinear meshes:
---------------------------

The `mesh::x(:,:)` array always contains a regular mesh, which gets 'distorded' to a curvilinear mesh by additional function calls.


Domain decomposition
--------------------

See [Domain_Decomposition.md](Domain_Decomposition.md).



