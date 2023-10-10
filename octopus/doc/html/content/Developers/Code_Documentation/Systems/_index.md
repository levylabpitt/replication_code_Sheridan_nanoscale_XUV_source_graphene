---
Title: "Systems in Octopus"
Weight: 11
---

{{< notice warning >}}
Work in progress!
{{< /notice >}}


The {{<developers "Code_Documentation:Multisystem_framework" "multisystem framework">}} allows to define and couple different physical systems.

The following diagram represents the family tree of the system classes. Rounded boxes denote abstract classes, while rectangular boxes are normal classes, which can be instantiated.

{{% graphviz-file "static/graph_data/interaction_partner_t.viz" %}}

### Abstract classes


#### {{< code interaction_partner_t >}}

An interaction_partner in {{< octopus >}} is anything, which can have an interaction with any other system.
For instance electrons and ions are interaction partners, but also the photonic system described by the Maxwell system, or external potentials.
Therefore, it is the base class of all possible systems and multisystems.

{{% expand "Definition of interaction_partner_t" %}}
```Fortran
#include_type_def interaction_partner_t
```
{{% /expand %}}

Each {{< name interaction_partner >}} is associated with a {{< name namespace >}}, owns a {{< name clock >}}, as well as a list of {{< developers "code_documentation/multisystems/interactions" "interactions" >}} in which it can be a partner, and a list of physical quantities, which can be exposed to other systems (through the interactions). See {{< developers "Code_Documentation:Multisystems:quantity" "here" >}} for the list of possible quantities.

It also provides the basic functions to update exposed quantities, and copy them to the interaction. More details about this mechanism are described in the section on 
{{< developers "code_documentation/multisystems/interactions" "interactions" >}}.

#### {{< code system_t >}}

The {{< code system_t >}} type is the abstract type for all systems. 
As all possible systems are potential partners of some interaction, the {{< code system_t >}} type itself extends the abstract {{< code interaction_partner_t >}} type.

{{% expand "Definition of system_t" %}}
```Fortran
#include_type_def system_t
```
{{% /expand %}}

The {{< code system_t >}} class adds information about the physical space in which the system exists, the propagator, and the list of interactions, which are owned by the system. (Check the section {{< developers "Code_Documentation:Multisystems/Interactions" "interactions" >}} for details on who owns an interaction.)


#### {{< code multisystem_t >}}

The {{< code multisystem_t >}}, finally adds a list of systems to the type definietion.

{{% expand "Definition of multisystem_t" %}}
```Fortran
#include_type_def multisystem_t
```
{{% /expand %}}

{{< notice note >}}
{{< code multisystem_t >}} is an abstract class and cannot be used as such to describe a set of systems in the code.
The {{< code type >}} to be used for combined systems is {{< code multisystem_basic_t >}}, which extends {{< code multisystem_t >}}.
{{< /notice >}}




### Specific classes

{{% children depth=2 %}}


Currently, the following system types are defined:
```Fortran
#include_code_doc system_types
```
{{% notice note %}}
When using these system types, **always** use the parameters, and not their numerical values, as they might change over time.
{{% /notice %}}




