---
Title: Mesh implementation
Weight: 2
---



The data structures:
--------------------

{{% expand "Definition of basis_set_abst_t" %}}
```Fortran
#include_type_def basis_set_abst_t
```
{{% /expand %}}

{{% graphviz-file "static/graph_data/basis_set_abst_t.viz" %}}


The mesh descriptor is defined in 
`grid/mesh.F90`:

```Fortran
#include_type_def mesh_t
```
which uses the structure `index_t`, containing the range and the mapping arrays:
```Fortran
#include_type_def index_t
```

About lxyz and lxyz_inv maps:
-----------------------------

The direct map:

> lxyz(_d_, _i_) = R<sup>_(d)_</sup><sub>_i_</sub>

where _d_ denotes the dimension (i.e. x, y or z) and _i_ is the index of the point.

The inverse map:

> lxyz_inv( R<sup>x</sup><sub>_i_</sub> , R<sup>y</sup><sub>_i_</sub> , R<sup>z</sup><sub>_i_</sub> ) 
= _i_

The points defined by lxyz define a rectangular box of _d_ dimensions. 
The real mesh vectors are related to R<sub>_i_</sub> by multiplication with spacing
and possibly distortion for curvilinear coordinates.

Note that this index array is the same in all parallel domains, and relates the integer coordinates to the global index of a given point.




```Fortran
#include_type_def mesh_cube_map_t
```







Note on packed states:
----------------------



Note on curvilinear meshes:
---------------------------

The `mesh::x(:,:)` array always contains a regular mesh, which gets 'distorded' to a curvilinear mesh by additional function calls.


Domain decomposition
--------------------



