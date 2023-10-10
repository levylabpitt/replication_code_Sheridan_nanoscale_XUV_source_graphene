---
Title: "Object orientation in Octopus"
Weight: 1
---


### Motivation 

Since version 10, {{< octopus >}} is gradually being converted into a fully object oriented code, using the OOP features of Fortran 2003.

The main features of OOP used in {{< octopus >}} are:

* Encapsulation
* Inheritance

The benefits of OOP design are:

* less code duplication
* better readable code
* code can be closer to the physics
* low level aspects can be hidden


### What is an object?

An object refers to a data-structure, which is bundled with the functions (methods) acting on that data. Usually, it is good practice to declare the data as private, 
and allow access only through so-called access functions (getters and setters). This is called encapsulation, and allows later to change details of the implementation
without affecting the rest of the code.

Objects are instances of classes, which define the data structures and the methods for its objects.


### Inheritance and class hierarchy 

Object oriented programming allows to build more complex classes by inheriting the structure and behaviour of other classes.
This reduces duplicated code, and ensures consistency between different classes, which share some functionality.

{{% expand "Example: linked list" %}}
Many classes require the functionality to iterate over a list of objects. These objects, however depend on the nature of the specific class.
In the diagram below, you see the inheritance tree of the {{<code linked_list_t >}} class, which is the parent for many more specialized classes.
{{% graphviz-file "static/graph_data/linked_list_t.viz" %}}
{{% /expand %}}



### Templating

Languages like C++ allow to use templates to write functions or algorithms independent of the data type, they are to be applied to.
Unfortunately, this is not directly possible in Fortran. {{< octopus >}} uses a poor-mans approach to templates, by using macros and including 'templates' of functions, where the explicit types are replaced by macros. These routines are usually defined in files which have '_inc' in their name.

To create a 'templated' function, one can either provide an explicit interface
```Fortran
interface my_function
  module procedure dmy_function, zmy_function
end interface
```
of the function needs to be called using the `X(...)` macro. 

The body of the function is usually defined in a separate file (here {{< code "my_funtion_inc.F90" >}}) which needs to be included as shown here:

```Fortran

#include "undef.F90"
#include "real.F90"
#include "my_function_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "my_function_inc.F90"
...
```

In {{< code my_function_inc.F90 >}} we would have something like:
```Fortran
function X(my_function)(arg1, arg2) result(res))

  R_TYPE, intent(in)  :: arg1
  R_TYPE, intent(in)  :: arg2
  R_TYPE, intent(out) :: res

  ...

end function X(my_function)
```

The used macros are defined in these include files:

{{% expand "include/real.F90" %}}
```Fortran
#include_file src/include/real.F90
``` 
{{% /expand %}}

{{% expand "include/complex.F90" %}}
```Fortran
#include_file src/include/complex.F90
``` 
{{% /expand %}}

{{% expand "include/integer.F90" %}}
```Fortran
#include_file src/include/integer.F90
``` 
{{% /expand %}}

{{% expand "include/undef.F90" %}}
```Fortran
#include_file src/include/undef.F90
``` 
{{% /expand %}}
