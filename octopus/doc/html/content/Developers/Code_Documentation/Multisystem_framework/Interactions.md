---
Title: Interactions
section: Developers
Weight: 2
---

{{< notice warning >}}
Work in progress!
{{< /notice >}}

### Introduction

In physics, interactions are in general symmetric. However, the magnitude with which it affects the interaction partners can be very different.
Also, the way the interactions have to be implemented can differ depending on the ''direction'' of the interaction.

Therefore, in {{< octopus >}} interactions are broken down into their forward and their backward parts, which are treated seperately.


{{% graphviz %}}
digraph {
  "partner A" -> "partner B" [label=" forward   "];
  "partner B" -> "partner A" [label=" backward   "];
 }
{{% /graphviz %}}

As we will see, in some approximations, only one of them (usually the forward interaction) are used.

The ''partners'' are instances of {{< code interaction_partner_t >}} derived types. The interactions (from now on uni-directional) are implemented as derived types of the abstract {{< code interaction_t >}} class.

{{% graphviz %}}
digraph {
  "partner A" -> "partner B" [label=" interaction "];
 }
{{% /graphviz %}}

An interaction ''belongs'' to partner A and acts on partner B.

### Abstract classes

#### {{< code interaction_t >}}

{{% expand "Definition of interaction_t" %}}
```Fortran
#include_type_def interaction_t
```
{{% /expand %}}

This is the minimal data structure for interactions, which only contains information relating to partner A, who owns the interacion.
In particular, the interaction has a list of quantities, which partner A needs to expose, so that the interaction can be calculated.
The IDs for the exposed quantities are defined in the section {{< developers "Code_Documentation:Multisystems:quantity" "Exposed Quantities" >}}.

Furthermore, this abstract type already contains the clock for the interaction, and defines the interfaces for the deferred {{< code "update()" >}} 
and {{< code "calculate()" >}} routines.

{{% graphviz-file "static/graph_data/interaction_t.viz" %}}


#### {{< code interaction_with_partner_t >}}

Curretly, all interactions are derived from {{< code interaction_with_partner_t >}}, which extends {{< code interaction_t >}} and adds the information about partner B, from here on only called ''partner''.

{{% expand "Definition of interaction_with_partner_t" %}}
```Fortran
#include_type_def interaction_with_partner_t
```
{{% /expand %}}

This type adds a pointer to the interaction partner and a list of IDs of the quantities, the partner exposes to the interaction.
Furthermore, at this level, we can implement the {{< code "update()" >}} procedure.


### Abstract children

#### {{< code force_interaction_t >}}

The next level of specialization of interaction, are all interactions which create a force on the partner. Here we only add the actual force vector, acting on the partner system.
{{% expand "Definition of force_interaction_t" %}}
```Fortran
#include_type_def force_interaction_t
```
{{% /expand %}}

#### {{< code potential_interaction_t >}}

The next level of specialization of interaction, are all interactions which create a potential for the partner. 

{{% expand "Definition of potential_interaction_t" %}}
```Fortran
#include_type_def potential_interaction_t
```
{{% /expand %}}


#### {{< code density_interaction_t >}}

The next level of specialization of interaction, are all interactions which create a potential for the partner. 

{{% expand "Definition of density_interaction_t" %}}
```Fortran
#include_type_def density_interaction_t
```
{{% /expand %}}

</br>

### Specific classes:

Specific interaction classes extend the abstract ones. The most important element they add to the abstract classes is the information about the quantities, required to calculate the interaction. In case of the system. owning the interaction (system A), it is sufficient to keep pointers to the data, stored in the system itself. Thr reason is that the interaction is always updated by the propagator of the system A. For the partner system (system B), however, the interaction keeps a copy of the exposed quantities. This allows the partner system to continue the propagation beyond the time for which the quantities are requested, which might happen if the two systems are using different time steps.




#### Ghost interaction

{{% expand "Definition of ghost_interaction_t" %}}
```Fortran
#include_type_def ghost_interaction_t
```
{{% /expand %}}

#### Gravity

{{% expand "Definition of gravity_t" %}}
```Fortran
#include_type_def gravity_t
```
{{% /expand %}}

#### Coulomb force

{{% expand "Definition of coulomb_force_t" %}}
```Fortran
#include_type_def coulomb_force_t
```
{{% /expand %}}

#### Lorentz force

{{% expand "Definition of lorentz_force_t" %}}
```Fortran
#include_type_def lorentz_force_t
```
{{% /expand %}}


### Interaction factory

Instances of {{< code "interaction_t" >}} or derived types are, like systems, generated using a {{< versioned-link "developers/code_documentation/multisystems/factories/#interaction-factories" "factory" >}}.


Currently, the following interaction types are defined:
```Fortran
#include_code_doc interaction_types
```
{{% notice note %}}
When using these system types, **always** use the parameters, and not their numerical values, as they might change over time.
{{% /notice %}}


{{% graphviz-file "/static/graph_data/interaction_graph.dot" %}}