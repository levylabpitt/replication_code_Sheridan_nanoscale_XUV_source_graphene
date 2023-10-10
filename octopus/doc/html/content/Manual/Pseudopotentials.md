---
title: "Pseudopotentials"
#series: "Manual"
weight: 100
description: " "
---


{{< octopus >}} is capable of using different formats for pseudopotentials (see {{< variable Species >}} for more information). We include in the distribution some LDA Troullier-Martins pseudopotentials for some common atoms, but it is likely that you will need other pseudopotentials.

#### Repositories on the web 

These are places on the web that include pseudopotential files in a variety of formats. Note that not all formats are read by {{octopus}} but many are.

* [NNIN Virtual Vault for Pseudopotentials](https://www.nnin.org/nnin_comp_psp_vault.html) Contains an extensive list of pseudo-potential databases and a search engine. It also has links to pseudo-potential related tools.
* [ABINIT](https://www.abinit.org) has a quite complete [set of pseudopotentials](https://www.abinit.org/atomic-data-files). At this moment only <tt>.fhi</tt> files are supported.
* [QuantumEspresso](https://www.quantum-espresso.org/)'s format is also accepted. So you can use the <i>norm-conserving</i> pseudopotentials from their [repository](https://www.quantum-espresso.org/pseudopotentials).

#### Pseudopotential generators 

In order to generate your own pseudopotentials, we advise the use of the following programs:
* [APE](https://www.tddft.org/programs/APE) 
* [opium](https://opium.sourceforge.net/) 
* [FHI98PP](https://www.fhi-berlin.mpg.de/th/fhi98md/fhi98PP/) 

