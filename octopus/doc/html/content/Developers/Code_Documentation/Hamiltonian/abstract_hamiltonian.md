---
Title: "Abstract Hamiltonian class"
Weight: 1
Description: "How to contribute to the code."
---

This abstract class does not contains any real data or computational terms, but rather defines the abstract interface, which is inherited by the specific Hamiltonian classes, and some information, which is common to all systems, such as whether the Hamiltonian is Hermitian, and some variables describing the spectral range.

```Fortran
#include_type_def hamiltonian_abst_t
```

This abstract class is then specialized for the various systems: 
{{% graphviz-file "static/graph_data/hamiltonian_abst_t.viz" %}}
