---
Title: "Linked list"
Weight: 1
---


The {{< code "linked_list_t" >}} type implements a linked list of unlimited polymorphic
values. This allows the storage of any type of data. 

{{% expand "Definition of linked_list_t" %}}
```Fortran
#include_type_def linked_list_t
```
{{% /expand %}}

where the data is stored in {{< code "list_node_t" >}} type objects.

{{% expand "Definition of list_node_t" %}}
```Fortran
#include_type_def list_node_t
```

{{% notice note %}}
The {{<code "class(*), pointer :: value => null()">}} allows storage of a pointer to {{<emph any>}} data type. 
{{% /notice %}}
{{% /expand %}}


Iterating over the list is done using the associated iterator. 
{{% expand "Definition of linked_list_iterator_t" %}}
```Fortran
#include_type_def linked_list_iterator_t
```
{{% /expand %}}


These classes are not meant to used as is, but rather to be extended and by providing an add method to the
list and a get_next method to the iterator.

{{% graphviz-file "static/graph_data/linked_list_t.viz" %}}

Also the {{< code "linked_list_iterator_t" >}} is the parent of a number of specialized iterators for the various systems:

{{% graphviz-file "static/graph_data/linked_list_iterator_t.viz" %}}
