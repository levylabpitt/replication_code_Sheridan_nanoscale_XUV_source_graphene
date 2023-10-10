---
title: "Updating to a new version"
---


This page lists the issues that may appear when updating to a new version of {{< octopus >}} and how to solve them.

####  Octopus 13  

* Libxc version 4 has been marked deprecated and will be removed in the next major version of Octopus.

#####  Variables 

* The variables TDIonicTimeScale and TDTimeStep are now disconected variables. TDTimeStep specifies the electronic time step, and the ionic timestep is given by TDTimeStep*TDIonicTimeScale

####  Octopus 11  

#####  Variables  

* The syntax of the Output and TDOutput variables has been changed for greater flexibility.
* The LSize variable cannot be used anymore for periodic systems. Instead, the LatticeParameters variable is mandatory for periodic systems.

#####  Execution  

* Runs that do not use the GPU will now by default fail if the code was compiled with GPU support. It is possible to override this behavior using the AllowCPUonly input variable.

#####  Installation  

* This version supports {{<libxc>}} 5.1. Support for {{<libxc>}} 3 and {{<libxc>}} 5.0 has been removed. Note that if you are using {{<libxc>}} 5.1, the library needs to be compiled with support for third order derivatives in order for Octopus to be able to perform response calculations using the Casida or the Sternheimer methods.

####  Octopus 10  

#####  Execution  

* The default convergence criteria have been changed: {{<code-inline >}} {{<variable "ConvRelDens">}} = 1e-6 {{</code-inline>}} and {{<code-inline>}}{{<variable "EigensolverTolerance">}} = 1e-7 {{</code-inline>}}.
* The SCF loop is terminated only if the convergence criteria are fulfilled twice in subsequent iterations.

####  Octopus 8  

#####  Variables  

* The {{< variable "Species" >}} block can now take a {{< emph "set" >}} option to choose which pseudopotential set to use.
* The default for {{< variable "XCFunctional" >}} is now taken, when possible, from the pseudopotentials.

####  Octopus 7  

#####  Variables  

* The UnitsInput and Units variables have been removed. Input file values are now always in atomic units, unless explicitly stated otherwise in the corresponding variable description. A few constants are available to help using eV/Angstrom units in the input file, such that one can write {{<code-inline>}} {{<variable "Spacing">}} = 0.25*angstrom {{</code-inline>}}.
* Now by default the code assumes XYZ files to be in Angstrom. This can be changed by using the {{< variable "UnitsXYZFiles" >}} variable.

####  Octopus 6  

#####  Installation  

* A compiler supporting Fortran 2003 iso_c_binding is now required.
* The executable is always called {{< file "octopus" >}}, no longer {{< file "octopus_mpi" >}} when compiled in parallel.
* This version supports {{<libxc>}} 3.0.0. Although it is still possible to use {{<libxc>}} 2.0.0 or any of the 2.x versions, updating to {{<libxc>}} 3.0.0 is highly recommended.

#####  Variables  

* There is a new format for the {{< variable "Species" >}} block and the names of the species options have been revised.
* The variable OutputHow has been renamed to {{< variable "OutputFormat" >}}.
* The variable ParallelizationStrategy and ParallelizationGroupRanks have been replaced with {{< variable "ParDomains" >}}, {{< variable "ParStates" >}}, {{< variable "ParKPoints" >}}, and {{< variable "ParOther" >}}.

#####  Utilities  

The {{< file "oct-rotatory_strength" >}} utility has been removed. This utility is replaced by a new {{< variable "PropagationSpectrumType" >}} in {{< file "oct-propagation_spectrum" >}}. 

####  Octopus 5.0  

#####  Restart files  

This version breaks the backwards compatibility of the restart files. Nevertheless, it is not too difficult to update the restart files generated with a previous version. Here is a summary of the changes and how to update the files:

* Casida restart data now goes to a directory {{< file "restart/casida" >}}, and a file {{< file "kernel" >}} or {{< file "kernel_triplet" >}}. One then needs to create the directory and rename the files.
* The format of the {{< file "occs" >}} file changed with the addition of an extra field for the imaginary part of the eigenvalues. To update the file one needs to add a column of zeros by hand after the column of eigenvalues.
* The line starting with {{< emph "fft_alpha" >}} was removed from the {{< file "mesh" >}} file.

Note that the mesh partition restart is also not compatible any more. In this case the partition needs to be recalculated, which is usually not a problem.

#####  Variables  

Some variables changed name or changed format. If an obsolete variable is found, {{< octopus >}} will stop and will tell you the variable that replaces it.

The datasets feature was removed. The input file of a calculation using datasets needs to be split into several independent input files.

####  Octopus 4.1  

#####  Installation  

This version requires {{<libxc>}} 2.0.0 or higher, so it might be necessary to update {{<libxc>}} before installing {{< octopus >}}. 

#####  Restart files  
Casida restart file is incompatible with that generated by version 4.0. The new format avoids problems with restarting with a different number of states.

####  Octopus 4.0  

#####  Installation  

* Libxc is now an independent library. To compile {{< octopus >}} 4.0.0 you will have to compile {{<libxc>}} 1.1.0 first. These are the short instructions to compile it (we assume that {{<libxc>}} will be installed in $DIR, this can be the same directory where {{< octopus >}} is going to be installed):
{{< code-block >}}
 tar -xvzf libxc-1.1.0.tar.gz
 cd libxc-1.1.0 
 ./configure --prefix=$DIR
 make
 make install
{{< /code-block >}}
Now, when configuring {{< octopus >}} pass the option --with-libxc-prefix=$DIR.

* The configure option for the location of netcdf and etsf_io are --with-netcdf-prefix and --with-etsf-io-prefix.

#####  Variables  

* TDEvolutionMethod is now called {{< variable "TDPropagator" >}}.

#####  Utilities  

* {{< file "oct-cross_section" >}} was renamed to {{< file "oct-propagation_spectrum" >}}.
* {{< file "oct-broad" >}} was renamed to {{< file "oct-casida_spectrum" >}}.
* The format for {{< file "oct-help" >}} has changed. Now {{< file "-s" >}} is used instead of {{< file "search" >}}, {{< file "-p" >}} instead of {{< file "show" >}}, and  {{< file "-l" >}} instead of {{< file "list" >}}.

####  Octopus 3.0  

#####  Input variables  

Some variables changed name or changed format. If an obsolete variable is found, {{< octopus >}} will stop and will tell you the variable that replaces it.

* {{< variable "Units" >}}, {{< variable "UnitsInput" >}} and {{< variable "UnitsOutput" >}} now take a named option as argument instead of a string. So

{{< code-block >}}
 Units = "eVA"
{{< /code-block >}}

should be replaced by
{{< code-block >}}
 
 Units = eV_Angstrom
{{< /code-block >}}

The code will stop if the old format is encountered.

* XFunctional and CFunctional were replaced by {{< variable "XCFunctional" >}}, for example

{{< code-block >}}
 XFunctional = lda_x
 CFunctional = lda_c_vwn
{{< /code-block >}}

must be replaced with

{{< code-block >}}
 XCFunctional = lda_x + lda_c_vwn
{{< /code-block >}}

* TDLasers was replaced by {{< variable "TDExternalFields" >}}.
* Some options for {{< variable "CalculationMode" >}} were renamed.

#####  Output directories  

Some directories in {{< octopus >}} output were renamed, restart files now are stored under {{< file "restart/" >}} instead of {{< file "tmp/" >}}. Files previously found under {{< file "status/" >}} are now located in {{< file "exec/" >}}.

#####  Restart file format  

{{< octopus >}} writes restart files in a binary format, this format has been updated and improved. As a result {{< octopus >}} is no longer able to restart automatically the calculation of a run performed with older versions. 

#####  Recovering old restart files  

If you really need to recover your old restart files, first you have to generate NetCDF restart files with your old version of {{< octopus >}}. Then you will have to rename the {{< file "tmp/" >}} directory to {{< file "restart/" >}} and remove the {{< file "restart_" >}} prefix from its subdirectories, for example {{< file "tmp/restart_gs/" >}} now should be called {{< file "restart/gs/" >}}.

Now you can run {{< octopus >}} 3.0 and it will find your restart information.


{{< manual-foot prev="Manual:Deprecated Utilities" next="Manual:Building from scratch" >}}
