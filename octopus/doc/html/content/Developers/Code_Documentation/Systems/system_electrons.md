---
Title: "System: Electrons"
Weight: 2
---

{{% notice warning %}}
The {{<code "electron_t">}} is currently in the process of being refactored. In the end, this class shall describe the electrons alone, which are interacting
with ions, external fields, etc. through interactions. Also the electron-electron interaction, described according to the various theory levels (see {{<variable "TheoryLevel">}}) are implemented through a so-called intra-interaction.
{{% /notice %}}

{{% expand "Definition of \"electrons_t\"" %}}
```Fortran
#include_type_def electrons_t
```
{{% /expand %}}

{{% expand "Definition of the constructor" %}}
```Fortran
#include_function electrons_constructor
```
{{% /expand %}}
