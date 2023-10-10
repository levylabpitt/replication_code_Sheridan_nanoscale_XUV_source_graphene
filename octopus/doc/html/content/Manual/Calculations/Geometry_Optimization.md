---
title: "Geometry Optimization"
#series: "Manual"
weight: 6
---


To perform a geometry optimization with {{< octopus >}} one should set the {{< variable "CalculationMode" >}} to {{< code "go" >}}. The method to be used should be set using the variable {{< variable "GOMethod" >}}. Note that most of the methods use the minimization routines from [GSL](https://www.gnu.org/software/gsl/manual/html_node/Multidimensional-Minimization.html).

The stopping criteria can be set with the variables {{< variable "GOTolerance" >}} and {{< variable "GOMinimumMove" >}}. Then minimization will be stopped when all forces on ions are smaller than {{< variable "GOTolerance" >}} or when all species coordinates change less than {{< variable "GOMinimumMove" >}} during one minimization step. If none of the previous criteria is matched after a number of minimization steps equal to {{< variable "GOMaxIter" >}}, then the minimization will be stopped with an error message.

After each minimization step taken by the [GSL](https://www.gnu.org/software/gsl/) algorithms {{< octopus >}} will write the current geometry to a file named {{< file "go.XXXX.xyz" >}} (XXXX is the iteration number) in the directory {{< file "geom" >}}. As an extra information, the title of the xyz file (second line) contains the total energy.

Currently, the default method used for geometry optimizations is {{< variable "GOMethod" >}} = {{< code "FIRE" >}}. Otherwise, many times, if the stopping criteria are too small, the GSL routines will return an error message stating that they were unable to improve the solution. 

At the end of the run, if the minimization was successful, the minimized geometry will be written to a file named {{< file "min.xyz" >}} in the working directory.

See also the tutorial on {{< tutorial "Unsorted/Geometry_optimization" "geometry optimization">}}

{{< manual-foot prev="Calculations:Optimal Control" next="Visualization" >}}
