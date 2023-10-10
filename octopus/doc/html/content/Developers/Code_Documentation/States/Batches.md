---
Title: Batches
Weight: 2
---


In many situations, we need to perform the same operations over many mesh functions, such as the electronic wave functions. 
It is therefore advantageous to group those functions into one object. This can ensure that different mesh functions are contiguous in memory.

Due to the nature of stencil operations, which constitute a large part of the low level operations on mesh functions, it is often more efficient to perform the same stencil operation over different mesh functions (i.e. using the state index as fast index), than looping first over the mesh index, which would, in general, require a different stencil for each mesh point. This is, in particular, the case for calculations utilizing GPUs.

Therefore, we store mesh functions in linear or in so-called packed form. The former refers to the 'natural' ordering where the mesh index is the fastest moving, while the latter is transposed. 

The abstract class {{< code "batch_t" >}} is the parent class for batches, such as electronic {{< versioned-link "Developers/Code_Documentation/States/States" "wave functions" >}}. 

{{% expand "Definition of \"batch_t\"" %}}
```Fortran
#include_type_def batch_t
```
{{% /expand %}}

This class includes information about the dimensions of the functions (number of states, spatial dimension and number of mesh points), but also internal book-keeping variables, 
keeping track of the status of the batch. Furthermore, the {{< code "batch_t" >}} data type contains pointers to the actual data arrays, and defines the methods for interacting with a batch.

Empty batches can be initialized with:
{{% expand "Initializing empty batches" %}}
```Fortran
#include_subroutine X(batch_init)
```
{{% /expand %}}



{{% expand "Initializing batches with memory" %}}
```Fortran
#include_subroutine X(batch_init_with_memory_1)
```

```Fortran
#include_subroutine X(batch_init_with_memory_2)
```

```Fortran
#include_subroutine X(batch_init_with_memory_3)
```
{{% /expand %}}


{{% graphviz-file "static/graph_data/batch_t.viz" %}}


{{% expand "Definition of \"wfs_elec_t\"" %}}
```Fortran
#include_type_def wfs_elec_t
```
{{% /expand %}}
