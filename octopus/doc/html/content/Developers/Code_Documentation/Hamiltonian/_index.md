---
Title: "Hamiltonian"
Weight: 20
Description: "How to contribute to the code."
---

{{< notice warning >}}
Work in progress!
{{< /notice >}}


The Grid in Octopus
====================

Introduction
------------


The fundamental element of Octopus is the real space grid, on which
all quantities, such as wafe functions, densities and potentials are
defined.

In general, a grid consists of

* the mesh
* the simulation box

The mesh is based on an regular mesh in _D_ dimensions, which can be distorted 
using so-called curvi-linear coordinates.

The underlying internal structure is a mesh of integer coordinates and the important
mapping arrays which 
* map the index to the _D_ integer coordinates
* map the integer coordinates to the index.

The grid is constructed in several steps:

1) Generation of a regular mesh filling a rectangular volume enclosing the simulation box.
2) Selecting inner points, enlargment points and ghost points. 
3) re-ordering the points and regeneration of the mapping arrays.



