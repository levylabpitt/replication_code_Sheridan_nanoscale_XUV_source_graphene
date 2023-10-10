---
title: "oct-convert"
#series: "Manual"
---


### Name 
oct-convert - Octopus utility to read obj files and write in many different file formats

### Description 

This executable gives the ability to read files written during the ground-state or time-dependent execution

### Example 

You can run ground-state and time-dependent execution of the [benzene example](../Benzene_molecule).

Then, we have to add this to the inp file, if we want to have the ground state density in DX format:

```bash
 {{< variable "Output" >}} = density
 {{< variable "OutputFormat" >}} = dx
 {{< variable "ConvertFolder" >}} = 'restart/gs'
 {{< variable "ConvertFilename" >}} = 'density'
 {{< variable "ConvertIterateFolder" >}} = no
```

To convert the restart wave-functions (from 1 to 10) of a td run:

```bash
 {{< variable "Output" >}} = density
 {{< variable "OutputFormat" >}} = dx
 {{< variable "ConvertIterateFolder" >}} = no
 {{< variable "ConvertFilename" >}} = ' '
 {{< variable "ConvertStart" >}} = 1
 {{< variable "ConvertEnd" >}}   = 10
 {{< variable "ConvertFolder" >}} = 'restart/td/'
```

If we want to convert the densities of the time-dependent executions, from files td.0000001 to td.0000010:

```bash
 {{< variable "Output" >}} = density
 {{< variable "OutputFormat" >}} = dx
 {{< variable "ConvertIterateFolder" >}} = yes
 {{< variable "ConvertStart" >}} = 1
 {{< variable "ConvertEnd" >}}   = 10
```


{{< manual-foot prev="Manual:External utilities:oct-conductivity" next="Manual:External utilities:oct-dielectric-function" >}}
---------------------------------------------
