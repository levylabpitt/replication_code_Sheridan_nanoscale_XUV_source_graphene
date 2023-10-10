---
title: "Profiling"
weight: 6
section: "Tutorial Profiling"
series: "Tutorials"
tutorials: ["HPC"]
difficulties: "Expert"
author: "Sebastian Ohlmann"
description: "Finding the hotspots of the code"
---


In this tutorial, we want to use the internal profiling capabilities of
{{< octopus >}} to check which parts of the code take most of the time for different
examples.

##  Profiling

Profiling means to measure certain performance characteristics for different
parts of a code. In {{< octopus >}}, we have implemented an internal profiling tool and its most important
capability is to measure the time spent in different regions of the code.
As it creates only negligible overhead, it is a good idea to always turn it on
to check if the simulation works as expected.

In the input file, the profiling can be controlled with the input variable
{{< variable "ProfilingMode" >}}. It has several options:
* '''prof_time''': Profile the time spent in defined profiling regions.
* '''prof_memory''': As well as the time, summary information on memory usage and the largest arrays are reported.
* '''prof_memory_full''': As well as the time and summary memory information, a log is reported of every allocation and deallocation. Unlike the other options, this one can create an non-negligible overhead.
* '''prof_io''': Count the number of file open and close.

The most important and useful option is usually {{<variable "ProfilingMode">}} {{<code "= prof_time">}} . 
The overhead for the memory profiling is bigger, so this should only be used when needed.

In the code, it is used as follows: first, an object of type {{<code "profile_t">}}
is created and then profiling_in/profiling_out is called. It looks like this:

```Fortran
use profiling_oct_m
...
subroutine ...
  type(profile_t), save :: prof
  call profiling_in(prof, "MY_FUNCTION")
  ...
  call profiling_out(prof)
end subroutine
```

There are few important points to note:
* The profile object (here ''prof'') must have the ''save'' argument. Otherwise, the result is not accumulated for all the calls of the profiled region.
* If you create such a region, be aware that there should not be two regions with the same name.

Using profiling is important when developing new code or when trying to optimize it, to keep track of the performance change and to know where to optimize.
The regions to be profiled are not necessarily entire routines and can correspond to only few lines of code, when this is relevant. 
The code can contains as many nested profiled region as needed, which are internally represented by a tree graph. This leads to the definition of self time and cumulative time, as we will discuss below.

##  Ground-state example

Let us first run a ground-state calculation with profiling. We will use
the following input file:

{{< code-block >}}
  {{< variable "CalculationMode" >}} = gs
  {{< variable "FromScratch" >}} = yes

  {{< variable "XYZCoordinates" >}} = "1ala.xyz"

  {{< variable "Radius" >}} = 4.0*angstrom
  {{< variable "Spacing" >}} = 0.4*angstrom

  {{< variable "ProfilingMode" >}} = prof_time
{{< /code-block >}}

In this example, the geometry of the molecule is specified in an extra file {{< file "1ala.xyz" >}}:

{{< code-block >}}
  23
 units: A
      N                   -1.801560    0.333315   -1.308298
      C                   -1.692266    1.069227    0.012602
      C                   -0.217974    1.151372    0.425809
      O                    0.256888    2.203152    0.823267
      C                   -2.459655    0.319513    1.077471
      H                   -1.269452    0.827043   -2.046396
      H                   -1.440148   -0.634968   -1.234255
      H                   -2.791116    0.267637   -1.602373
      H                   -2.104621    2.114111   -0.129280
      H                   -2.391340    0.844513    2.046396
      H                   -2.090378   -0.708889    1.234538
      H                   -3.530691    0.246022    0.830204
      N                    0.476130   -0.012872    0.356408
      C                    1.893957   -0.046600    0.735408
      C                    2.681281    0.990593   -0.107455
      O                    3.486946    1.702127    0.516523
      O                    2.498931    1.021922   -1.333241
      C                    2.474208   -1.425485    0.459844
      H                    0.072921   -0.880981    0.005916
      H                    1.975132    0.211691    1.824463
      H                    1.936591   -2.203152    1.019733
      H                    3.530691   -1.461320    0.761975
      H                    2.422706   -1.683153   -0.610313
{{< /code-block >}}

Now run the code. Alongside with the usual {{< octopus >}} folders, you will have a directory called `profiling`, that contains the files {{<file "time.000000">}} and {<<file "time.000000.tree<">}}.

The first file should look similar to the following output:

{{< code-block >}}
                                                                    CUMULATIVE TIME                                 |                         SELF TIME
                                          --------------------------------------------------------------------------|-------------------------------------------------------------
TAG                           NUM_CALLS      TOTAL_TIME   TIME_PER_CALL        MIN_TIME    MFLOPS  MBYTES/S   %TIME |       TOTAL_TIME   TIME_PER_CALL    MFLOPS  MBYTES/S   %TIME
====================================================================================================================|=============================================================
dNL_OPERATOR_BATCH                15544        2.721022        0.000175        0.000059    3671.2       0.0    28.9 |         2.721022        0.000175    3671.2       0.0    28.9
PS_FILTER                             4        1.552428        0.388107        0.266518       0.0       0.0    16.5 |         1.552428        0.388107       0.0       0.0    16.5
dDOTPV_MF_BATCH                   40017        1.319679        0.000033        0.000018    3918.5       0.0    14.0 |         1.319679        0.000033    3918.5       0.0    14.0
dBATCH_AXPY_FUNCTION              40017        1.165627        0.000029        0.000012       0.0       0.0    12.4 |         1.165627        0.000029       0.0       0.0    12.4
dVNLPSI_MAT_BRA                    8296        0.320456        0.000039        0.000034     827.2       0.0     3.4 |         0.320456        0.000039     827.2       0.0     3.4
dMF_DOTP                          54182        0.282613        0.000005        0.000003    6728.5       0.0     3.0 |         0.282613        0.000005    6728.5       0.0     3.0
EIGEN_SOLVER                         21        7.129506        0.339500        0.146633    2951.6     328.1    75.8 |         0.210328        0.010016    5802.8       0.0     2.2
dGET_POINTS                        1512        0.199078        0.000132        0.000021       0.0       0.0     2.1 |         0.199078        0.000132       0.0       0.0     2.1
dVLPSI                             8296        0.149986        0.000018        0.000013    2182.7   15731.6     1.6 |         0.149986        0.000018    2182.7   15731.6     1.6
dVNLPSI_MAT_KET                    8296        0.244048        0.000029        0.000026     670.0       0.0     2.6 |         0.129980        0.000016    1078.2       0.0     1.4
dPROJ_MAT_SCATTER                 91256        0.114068        0.000001        0.000001     204.8       0.0     1.2 |         0.114068        0.000001     204.8       0.0     1.2


...
{{< /code-block >}}

Each line corresponds to one profiling region in the code. The left half of the
table gives the cumulative time of each region (containing the time of all
sub-regions) and the number of calls to the region and the right half of the
table shows the self-time for each region. The rows are ordered by the
self-time. There are also some entries on MFLOPS and MBYTES/S, but those are
from analytical formulae in the code and they are not available for all regions
and they are also not always correct.

For this example, you can see that most of the time is spent in
"dNL_OPERATOR_BATCH" (29%), which is the non-local operator, i.e. the
application of the finite-difference stencil for computing derivatives like the Laplacian. This is good because it is usually
expected to be the most expensive part of the code. The "d" prefix means that this applied to real wavefunctions, whereas "z" would indicate that this applies to complex numbers. The other lines correspond
to other parts of the code.

Now change the eigensolver to RMMDIIS. For this, you need to set  {{<variable "Eigensolver" >}}{{<code  " = rmmdiis">}} and also set the number of {{<variable "ExtraStates" >}} to at least 10-20% of the total
number of states, let's choose 6 here (we have 32 states). So add {{<variable "ExtraStates" >}}{{<code " = 6">}} to the input file and run the code again. The output should be similar to:

{{< code-block >}}

                                                                    CUMULATIVE TIME                                 |                         SELF TIME
                                          --------------------------------------------------------------------------|-------------------------------------------------------------
TAG                           NUM_CALLS      TOTAL_TIME   TIME_PER_CALL        MIN_TIME    MFLOPS  MBYTES/S   %TIME |       TOTAL_TIME   TIME_PER_CALL    MFLOPS  MBYTES/S   %TIME
====================================================================================================================|=============================================================
dDOTPV_BATCH                      15320        2.005789        0.000131        0.000082    1018.6       0.0    21.6 |         2.005789        0.000131    1018.6       0.0    21.6
PS_FILTER                             4        1.947803        0.486951        0.338513       0.0       0.0    20.9 |         1.947803        0.486951       0.0       0.0    20.9
dNL_OPERATOR_BATCH                 4816        1.651675        0.000343        0.000097    7710.2       0.0    17.8 |         1.651675        0.000343    7710.2       0.0    17.8
dGET_POINTS                        4060        0.647682        0.000160        0.000047       0.0       0.0     7.0 |         0.647682        0.000160       0.0       0.0     7.0
dBATCH_AXPY_VEC                    5290        0.375005        0.000071        0.000046    1881.3       0.0     4.0 |         0.375005        0.000071    1881.3       0.0     4.0
dVLPSI                             3486        0.269985        0.000077        0.000017    1701.6    8219.7     2.9 |         0.269985        0.000077    1701.6    8219.7     2.9
dSET_POINTS                        1740        0.246499        0.000142        0.000058       0.0       0.0     2.7 |         0.246499        0.000142       0.0       0.0     2.7
dSET_BC                            4816        0.211553        0.000044        0.000007       0.0       0.0     2.3 |         0.211553        0.000044       0.0       0.0     2.3
dVNLPSI_MAT_BRA                    3486        0.204821        0.000059        0.000043    1218.0       0.0     2.2 |         0.204821        0.000059    1218.0       0.0     2.2
MESH_INIT                             3        0.120044        0.040015        0.000001       0.0       0.0     1.3 |         0.120044        0.040015       0.0       0.0     1.3
SG_PCONV                             31        0.114705        0.003700        0.003605       0.0       0.0     1.2 |         0.114705        0.003700       0.0       0.0     1.2
dSUBSPACE_HAMILTONIAN                29        0.574413        0.019807        0.019530    4312.0     324.4     6.2 |         0.108895        0.003755   13141.2       0.0     1.2
dSTATES_ROTATE                       29        0.366512        0.012638        0.012494    3904.4       0.0     3.9 |         0.107754        0.003716   13280.3       0.0     1.2
dRMMDIIS                             24        4.153825        0.173076        0.170564    2977.5     371.3    44.7 |         0.086432        0.003601       0.0       0.0     0.9
...
{{< /code-block >}}

Here, you can see that dNL_OPERATOR_BATCH is called less often and needs less time
in total - which is due to the fact that RMMDIIS uses batch operations as opposed
to the default CG eigensolver. Most of the time (22%) is spent in dDOTPV_BATCH
which is used in RMMDIIS itself and in the subspace diagonalization.

Now run the same simulation with a smaller spacing of 0.2. How does the
profiling change for CG and for RMMDIIS? Which operations would you expect
to take longer and which would you expect to be the same?

Now run the same simulation on 2 and on 4 cores. How does the profiling change?

##  TD example

Use the GS from the previous section as a starting point for analyzing TD run.
For this, run the input file

{{< code-block >}}
  {{< variable "CalculationMode" >}} = td
  {{< variable "FromScratch" >}} = yes

  {{< variable "XYZCoordinates" >}} = "1ala.xyz"

  {{< variable "Radius" >}} = 4.0*angstrom
  {{< variable "Spacing" >}} = 0.4*angstrom

  {{< variable "ProfilingMode" >}} = prof_time

  {{< variable "TDPropagator" >}} = aetrs
  {{< variable "TDMaxSteps" >}} = 20
  {{< variable "TDTimeStep" >}} = 0.05
{{< /code-block >}}

Then run the code. The {{<file "profiling/time.000000">}} file should be similar to

{{< code-block >}}
                                                                    CUMULATIVE TIME                                 |                         SELF TIME
                                          --------------------------------------------------------------------------|-------------------------------------------------------------
TAG                           NUM_CALLS      TOTAL_TIME   TIME_PER_CALL        MIN_TIME    MFLOPS  MBYTES/S   %TIME |       TOTAL_TIME   TIME_PER_CALL    MFLOPS  MBYTES/S   %TIME
====================================================================================================================|=============================================================
PS_FILTER                             4        1.600412        0.400103        0.267230       0.0       0.0    47.9 |         1.600412        0.400103       0.0       0.0    47.9
zNL_OPERATOR_BATCH                 1630        0.733867        0.000450        0.000286   14810.9       0.0    21.9 |         0.733867        0.000450   14810.9       0.0    21.9
zVLPSI                             1630        0.233238        0.000143        0.000114    3728.1    8046.5     7.0 |         0.233238        0.000143    3728.1    8046.5     7.0
MESH_INIT                             3        0.108109        0.036036        0.000001       0.0       0.0     3.2 |         0.108109        0.036036       0.0       0.0     3.2
zVNLPSI_MAT_BRA                    1630        0.095764        0.000059        0.000041    2968.5       0.0     2.9 |         0.095764        0.000059    2968.5       0.0     2.9
BLAS_AXPY_4                        1600        0.072772        0.000045        0.000021    8796.7       0.0     2.2 |         0.072772        0.000045    8796.7       0.0     2.2
zSET_BC                            1630        0.071694        0.000044        0.000023       0.0       0.0     2.1 |         0.071694        0.000044       0.0       0.0     2.1
SG_PCONV                             21        0.061658        0.002936        0.002843       0.0       0.0     1.8 |         0.061658        0.002936       0.0       0.0     1.8
BATCH_COPY_DATA_TO                 1200        0.058328        0.000049        0.000024       0.0       0.0     1.7 |         0.058328        0.000049       0.0       0.0     1.7
LIBXC                                84        0.042868        0.000510        0.000286       0.0       0.0     1.3 |         0.042868        0.000510       0.0       0.0     1.3
zPROJ_MAT_SCATTER                 17930        0.041486        0.000002        0.000002     747.7       0.0     1.2 |         0.041486        0.000002     747.7       0.0     1.2
COMPLETE_RUN                          1        3.344058        3.344058        3.344058    4017.2     569.1   100.0 |         0.039513        0.039513       0.0       0.0     1.2
...
{{< /code-block >}}

As you can see, we actually spent most of the time in PS_FILTER, which is the
filtering of the pseudopotentials. However, this is called only once in the
beginning for each species and thus its fraction reduces when running more
timesteps. Run the simulation for 200 timesteps and you will notice that its
share of the total time reduces. The next row shows zNL_OPERATOR_BATCH, which
is again the stencil operation, but this time for complex values (z...).

Now change the {{<variable "TDPropagator" >}} to exp_mid and run the simulation again (you
probably should not use it in practice). You will see in the profiling output
that the non-linear operator is only called about half as often as for aetrs,
as expected from the formula. The downside is of course that it's stability is
worse than aetrs, so usually you need larger timesteps for a stable propagation.

Next, try a different exponential method. Use {{<variable "TDExponentialMethod" >}}{{<code "= lanczos">}}
and set {{<variable "TDExpOrder" >}}{{<code "= 16">}} and run the simulation again. Do you notice
differences in the profiling? This exponential method can often be used for larger
timesteps than the default Taylor expansion, while still being numerically stable over long time propagations.


{{< tutorial-footer >}}


