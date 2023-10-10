---
Title: Factories
Weight: 6
---

Octopus uses the so-called [factory pattern](https://en.wikipedia.org/wiki/Factory_method_pattern) to create instances of the systems and interaction classes.
The abstract factory classes are introduced to avoid the problem of circular dependencies.

The function of the factories is to create an object of a (dynamically) given type and return a pointer to it. This is done by calling the respective 
constructors of the classes, of which an instance is to be created.

### System factories

```Fortran
#include_type_def system_factory_abst_t
```

{{% graphviz-file "static/graph_data/system_factory_abst_t.viz" %}}

```Fortran
#include_type_def system_factory_t
```

{{% expand "Definition of system_factory_create()" %}}
```Fortran
#include_function system_factory_create
```
{{% /expand %}}





### Interaction factories


```Fortran
#include_type_def interactions_factory_abst_t
```

{{% graphviz-file "static/graph_data/interactions_factory_abst_t.viz" %}}

```Fortran
#include_type_def interactions_factory_t
```

{{% expand "Definition of interactions_factory_abst_create_interactions()" %}}
```Fortran
#include_subroutine interactions_factory_abst_create_interactions
```
{{% /expand %}}


{{% expand "Definition of interactions_factory_create()" %}}
```Fortran
#include_function interactions_factory_create
```
{{% /expand %}}

