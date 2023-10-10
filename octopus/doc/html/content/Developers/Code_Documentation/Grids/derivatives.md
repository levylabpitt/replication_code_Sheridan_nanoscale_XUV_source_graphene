---
Title: Derivatives
Weight: 10
---

{{< octopus >}} use finite differences to evaluate derivatives, such as the Laplacian (kinetic energy).
The derivative at a point of the mesh is a weighted sum over neighboring points. For instance, 
the general form for the Laplacian is:
$$
    \nabla^2f(n_xh,\,n_yh) = \sum_i^n\sum_j^n\frac{c_{ij}}h\,f(n_xh + ih,\,n_yh+jh)
$$
The coefficients $\(c_{ij}\)$ depend on the mesh and number of points used and define the {{< emph stencil >}}.

The derivatives object contains operators for the gradient and the Laplacian, as well as information about the order and the stencil type, and information about boundaries:
```Fortran
#include_type_def derivatives_t
```

The derivative operators themselves are represented as non-local operators:
```Fortran
#include_type_def nl_operator_t
```

```Fortran
#include_type_def stencil_t
```

```Fortran
#include_type_def stargeneral_arms_t
```
