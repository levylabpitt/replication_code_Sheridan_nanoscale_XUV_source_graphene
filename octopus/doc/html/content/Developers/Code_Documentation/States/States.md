---
Title: Wave functions
Weight: 1
---

{{< notice warning >}}
Work in progress!
{{< /notice >}}

Wavefunctions in Octopus
========================

The wave functions in Octopus are referred to as the states.


They are handled by the module `states_abst_oct_m`, which is defined in [states/states_abst.F90](https://gitlab.com/octopus-code/octopus/blob/develop/src/states/states_abst.F90).

The exact way how states are stored in memory is flexible and depends on optimization  (and accelerator) settings in the input file. In particular, the order of indices depends on the `PACKED` setting.

States are stored in a hierarchy of 'containers'. Concepts of this hierarchy 
include _groups_, _batches_ and _blocks_.

#### The abstract class

The top level data structure, describing states is:

```Fortran
#include_type_def states_abst_t
```

This structure contains mainly metadata about states, describing _how_ states are represented, and defines the interface to the class.
As it is an abstract class, it cannot contain any information about the actual system. 

#### The states class for electrons

The class {{< code "states_elec_t" >}}, specialized the abstract class and contains more data specific to the electron system,
as well as pointers to other quantities, which are common to all states, such as the density, the current, etc.

{{% expand "Definition of \"states_elec_t\"" %}}
```Fortran
#include_type_def states_elec_t
```
{{% /expand %}}



The dimensions object contains a number of variables. The most relevant for this discussion is the `dim` variable, which denotes the dimension of one state, being `1` for spin-less states and `2` for spinors.

```Fortran
#include_type_def states_elec_dim_t
```


The wave functions themselves are stored in
```Fortran
    type(states_elec_group_t)     :: group
```
which, in turn, is defined in the module `states_elec_group_oct_m` in [`src/states/states_elec_group.F90`](https://gitlab.com/octopus-code/octopus/blob/develop/src/states/states_elec_group.F90):


The `group` contains all wave functions, grouped together in blocks or batches.
They are organised in an array of `batch_t` structures.

```Fortran
#include_type_def states_elec_group_t
```

```Fortran
    type(wfs_elec_t), pointer   :: psib(:, :)            !< A set of wave-functions
```
The indexing is as follows: `psib(ib,iqb)` where `ib` is the block index, and `iqn` the **k**-point. See below for the routine `states_init_block(st, mesh, verbose)` which creates the `group` object. On a given node, only wave functions of local blocks are available.
The `group` object does contain all information on how the batches are distributed over nodes.


```Fortran
#include_type_def wfs_elec_t
```




Creating the wave functions
---------------------------
A number of steps in initializing the `states_t` object are called from the `system_init()` routine:

`states_init()`:  
parses states-related input variables, and allocates memory for some book keeping variables. It does not allocate any memory for the states themselves.

`states_distribute_nodes()`:  
...


`states_density_init()`:  
allocates memory for the density (`rho`) and the core density (`rho_core`).

`states_exec_init()`:  
1. Fills in the block size (`st\%d\%block_size`);
2. Finds out whether or not to pack the states (`st\%d\%pack_states`);
3. Finds out the orthogonalization method (`st\%d\%orth_method`).
  


Memory for the actual wave functions is allocated in `states_elec_allocate_wfns()` which is called from the corresponding `*_run()` routines, such as `scf_run()` or `td_run()`, etc.


```Fortran
#include_subroutine states_elec_allocate_wfns
```

The routine `states_init_block` initializes the data components in `st` that describe how the states are distributed in blocks:

`st%nblocks`: this is the number of blocks in which the states are divided. 
Note that this number is the total number of blocks,
regardless of how many are actually stored in each node.  
`block_start`: in each node, the index of the first block.  
`block_end`: in each node, the index of the last block.
If the states are not parallelized, then `block_start` is 1 and `block_end` is `st%nblocks`.  
`st%iblock(1:st%nst, 1:st%d%nik)`: it points, for each state, to the block that contains it.  
`st%block_is_local()`: `st%block_is_local(ib)` is `.true.` if block `ib` is stored in the running node.  
`st%block_range(1:st%nblocks, 1:2)`: Block ib contains states fromn st\%block_range(ib, 1) to st\%block_range(ib, 2)  
`st%block_size(1:st%nblocks)`: Block ib contains a number st\%block_size(ib) of states.  
`st%block_initialized`: it should be .false. on entry, and .true. after exiting this routine.  

The set of batches `st%psib(1:st%nblocks)` contains the `block`s themselves.
  
```Fortran
#include_subroutine states_elec_init_block
```



The allocation of memory for the actual wave functions is performed in `batch_init_empty()` and `X(batch_allocate)()`. This routine, and the related `X(batch_add_state)()` show most clearly how the different memory blocks are related.

[`batch_init_empty()`](https://gitlab.com/octopus-code/octopus/blob/develop/src/grid/batch.F90) allocates the memory for `batch_state_t` `states` and `batch_states_l_t` `states_linear` and _nullifies_ the pointers within this types. Note that no memory for the actual wave functions has been allocated yet.

```Fortran
#include_subroutine batch_init_empty
```




Questions:
----------

How are the different objects pointing to states related?

* The usual storage for states is in `states_linear` which can be shadowed by `pack` in case 
  packed states are used.

What is the difference between `batch_add_state` and `batch_add_state_linear`?
