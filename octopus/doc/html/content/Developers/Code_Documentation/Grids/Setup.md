---
Title: Mesh Setup
Weight: 3
---


Setting up a grid
==================




In particular, we have:


* the "first stage":

  This sets up derivatives, double-grid and calls stages 1 and 2 of mesh_init.
At this stage we are still unaware or parallelism and domain decomposition.

{{% expand %}}
```Fortran
#include_subroutine grid_init_stage_1
```
{{% /expand %}}

* the second stage:

{{% expand %}}
```Fortran
#include_subroutine grid_init_stage_2
```
{{% /expand %}}


Mesh initialization:
--------------------

* First stage: set up the mesh_t::idx index structure, defining the lower and upper bounds.

{{% expand %}}
```Fortran
#include_subroutine mesh_init_stage_1
```
{{% /expand %}}


* Second stage: 

{{% expand %}}
```Fortran
#include_subroutine mesh_init_stage_2
```
{{% /expand %}}



* Third stage:

  We can finally generate the list of mesh points mesh_t::x.
At this stage, domain decomposition will be considered and x(:,:) contains
only the domain local points.

{{% expand %}}
```Fortran
! ---------------------------------------------------------
!> When running parallel in domains, stencil and np_stencil
!! are needed to compute the ghost points.
!! mpi_grp is the communicator group that will be used for
!! this mesh.
! ---------------------------------------------------------

#include_subroutine mesh_init_stage_3
```
{{% /expand %}}

which 'serializes' the mesh by a choice of space-filling curve methods.


