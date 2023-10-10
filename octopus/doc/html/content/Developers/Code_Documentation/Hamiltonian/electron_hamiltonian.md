---
Title: "Electron Hamiltonian class"
Weight: 2
Description: "How to contribute to the code."
---

The electronic Hamiltonian is derived from the abstract Hamiltonian. 

{{% expand "Definition of \"hamiltonian_elec_t\"" %}}
```Fortran
#include_type_def hamiltonian_elec_t
```
{{% /expand %}}

It contains the 'physical' quantities, such as 
- information about the dimensionality,
- the contributions to the potential, 
- a pointer to the ions,
- electronic mass,
- etc.

and an instance of {{< code "hamiltonian_elec_base_t" >}}, which bundles some lower level variables.
The separation of quantities into these two classes is mostly historic, and currently has no deeper systematics.

{{% expand "Definition of \"hamiltonian_elec_base_t\"" %}}
```Fortran
#include_type_def hamiltonian_elec_base_t
```
{{% /expand %}}

The method for application of the Hamiltonian, is implemented as a wrapper routine, which checks the compatibility of the supplied wave functions,
and then calls the batched version of the routine.
{{% expand "Applying the Hamiltonian (wrapper)" %}}
```Fortran
#include_subroutine X(hamiltonian_elec_apply)
```
{{% /expand %}}

This batched routine now takes care of the actual application of the various terms of the Hamiltonian. The different terms to be applied
can be selected using the optional {{< code "terms" >}} argument.
{{% expand "Applying the Hamiltonian (batched)" %}}
```Fortran
#include_subroutine X(hamiltonian_elec_apply_batch)
```
{{% /expand %}}

{{% notice "note" %}}
The Hamiltonian is only allowed to modify the wave functions. The reason why also the initial state {{< code "psib" >}} has {{< code "intent(INOUT)" >}} is that the Hamiltonian can pack the wave functions, if requested. 
{{% /notice %}}