---
title: "Scaling"
weight: 4
section: "Tutorial Scaling"
series: "Tutorials"
tutorials: ["HPC"]
author: "Sebastian Ohlmann"
description: "How many processors can I use?"
---


Octopus can be used on large clusters and supercomputers, but to do this
efficiently, it is important to understand how to judge the efficiency of its
parallelization scheme. To achieve this, we will conduct a
'''scaling analysis''' in this tutorial, for which we will run Octopus for a
certain system with different numbers of processes and compare the timings. This
kind of analysis is also called "strong scaling".


##  Introduction

A code like octopus is parallelized to utilize more computing resources at the
same time, thus reducing the time to solution for calculations. Naively, one
could expect that using twice the number of CPU cores for a calculation should
cut the time needed by a factor of two. However, not all parts of a program are
typically running in parallel and also there is some overhead associated with
the parallelization scheme itself, such as communication between processes.
Hence, the time $T(N)$ for the execution of a certain computation on
$N$ processors will depend on the fraction of serial work
$s$ because only the parallel fraction $1-s$ can be sped
up by a factor of $N$:

$$
T(N) = s T(1) + \\frac{1-s}{N}T(1)
$$

The ''speed-up'' $S$ is defined as the ratio $T(1)/T(N)$.
Plugging in the formula from above yields

$$
S = \\frac{T(1)}{T(N)} = \\frac{1}{s + \\frac{1-s}{N}} \\quad
\\xrightarrow{N\\to\\infty} \\frac{1}{s}
$$

This relation is called ''Amdahl's law''. It means that the achievable speed-up
can never be larger than the inverse of the serial fraction. If the serial
fraction is 50%, the speed-up cannot be larger than 2; if the serial fraction is
10%, the speed-up cannot be larger than 10; if the serial fraction is 1%, the
speed-up cannot be larger than 100. As you can see, the serial fraction must be
very small to achieve large speed-ups: for achieving a speed-up of 10000, the
serial fraction must be less than $10^{-4}$.

From this, it also follows that it one should not use more cores for a
computation than the maximum speed-up that can be reached. If, for example, the
maximum speed-up is 10, it does not make sense to run this calculation on 100
cores - this would be inefficient and waste resources.

To judge the efficiency of the parallelization, one can use the
''parallel efficiency'' which is defined as the ratio of the observed speed-up
to the ideal speed-up, $\epsilon = \frac{S}{S_\text{ideal}}$. As an
example, when comparing a run on 4 cores to a run on 1 core, for which the
speed-up is 3, the ideal speed-up would be 4 and thus the parallel efficiency is
$\epsilon = 3/4 = 75%$ in this case. As a rule of thumb, efficiencies
above 70% are acceptable.

In practice, one usually does not know the serial fraction. Thus, one executes
the code for a short test case for several numbers of processes to obtain the
speed-up, scaling curve, and efficiency. From this, one can then infer up to
which point the parallelization is still efficient and one should choose this
number of processors for subsequent runs to efficiently use the computing
resources.

When in doubt, it is usually more efficient to use slightly less resources, but
run several simulations at the same time (which is often required) - this will
use the resources efficiently and still provide a small total time to solution
because several calculations can run in parallel.

Be aware that on supercomputers, the available computing time is shared between
the members of a project or between the members of a group. If everyone runs
their code efficiently, the number of simulations that can be done by all
members of the group is maximized and thus more papers can be produced!


##  Scaling analysis for domain parallelization

As a first example, we will do a strong scaling analysis for the same input as
in the {{<tutorial "HPC/Parallelization" "previous tutorial">}}.
Save the following as {{< file "inp" >}} file:
```
#include_input doc/tutorials/hpc/scaling/1.domain_parallelization/inp
```
and the following as "1ala.xyz":
```
#include_file doc/tutorials/hpc/scaling/1.domain_parallelization/1ala.xyz
```

Then, run the ground state calculation or reuse the ground state from the
previous tutorial. For running the calculations, you can use the batch scripts
from the previous tutorial.

For this analysis, we will only look at the parallelization using MPI and
disregard the OpenMP parallelization for the time being. Usually, these two
components can be examined independent of each other, keeping the other fixed.

###  Baseline: serial run

As a baseline, we need to run the TD calculation on one processor. Change the
CalculationMode variable to td and submit a batch script to execute it on one
core. Check the parallelization section in the output and make sure that it is
executed in serial.

From the output of the timestep information, you can again get an average of the
elapsed time per timestep:

{{< code-block >}}
#include_file doc/tutorials/hpc/scaling/1.domain_parallelization/time-dependent.txt
{{< /code-block >}}

In this case, the average of the 20 time steps would be about 0.90s.


###  Scaling runs

To analyze the scaling, we will run TD calculations with a logarithmic spacing
in the number of processors, using 2, 4, 8, and 16 cores. To make sure that
domain parallelization is used, add the following to the {{< file "inp" >}} file:

{{< code-block >}}
ParStates = 1
ParDomains = auto
{{< /code-block >}}

Submit a job script with for each of the core numbers and wait until they are
finished. Then, extract the average time for one time step as we did above for
the serial run.

With these numbers, you can fill the following table:

|Cores          |  1   |  2   |  4   |  8   |  16  |
| :-----------: | :--: | :--: | :--: | :--: | :--: |
|Time           | 0.90 | 0.57 | 0.34 | 0.21 | 0.11 |

The numbers you get from your own measurement can deviate, but overall the
result should be similar.


###  Analysis

Using the times you measured, compute the speed-up and the parallel efficiency
as introduced earlier in the tutorial. You should get a table that is similar to
the following:


|Cores          |  1   |  2   |  4   |  8   |  16  |
| :-----------: | :--: | :--: | :--: | :--: | :--: |
|Time           | 0.90 | 0.57 | 0.34 | 0.21 | 0.11 |
|Speed-up       | 1.00 | 1.58 | 2.65 | 4.29 | 8.18 |
|Ideal Speed-up | 1    | 2    | 4    | 8    | 16   |
|Efficiency     | 1    | 0.79 | 0.66 | 0.54 | 0.51 |

What do these numbers tell us? First of all, the parallel efficiency is almost
80% when using 2 cores, so that is still fine. Going to more cores, the
efficiency drops. Thus, for this system it is inefficient to use more than 2
cores for the domain parallelization.

A standard way of visualizing scaling data is to display a log-log plot of the
speed-up vs. the number of cores. A little python script for doing this can be
found at the end of the tutorial. When the speed-up is near the line of the
ideal speed-up, the efficiency is good, but once it drops below, the scaling
breaks down.

Why does the efficiency drop for more than 2 cores? As you might remember from the
{{<tutorial "HPC/parallelization" "previous tutorial on parallelization" >}},
the parallelization in domains is only efficient, when the number of inner
points is large enough compared to the number of ghost points. So let's compare
the ratio of ghost points to local points, which should be less than about 25%
as a rule of thumb. For this, you need to look at the information on the mesh
partitioning, which should look similar to the following on 2 cores:

{{< code-block >}}
#include_file doc/tutorials/hpc/scaling/1.domain_parallelization/partition.txt
{{< /code-block >}}

For this case, the ratio would be about 14435/70175=21%.

What are the ratios you get for the other core numbers? Can they explain the
drop in efficiency?

(as a reference, the ratios should be roughly:
2 cores: 21%; 4: 35%; 8: up to 70% 16: up to 111%)


##  Scaling analysis for states parallelization

Let's do a similar analysis, but for states parallelization. Change the
parallelization options to

{{< code-block >}}
ParStates = auto
ParDomains = 1
{{< /code-block >}}

to make sure only states parallelization is used. Then submit batch jobs to run
the code on 1, 2, 4, 8, and 16 cores. Gather the timings as in the previous
section and create a table of the timings; also compute the speed-up and
parallel efficiency for all runs. The table should be similar to (the timings
can vary, but the speed-up and efficiency should be similar):

|Cores          |  1   |  2   |  4   |  8   |  16  |
| :-----------: | :--: | :--: | :--: | :--: | :--: |
|Time           | 0.67 | 0.36 | 0.23 | 0.14 | 0.15 |
|Speed-up       | 1.00 | 1.86 | 2.91 | 4.79 | 4.47 |
|Ideal Speed-up | 1    | 2    | 4    | 8    | 16   |
|Efficiency     | 1    | 0.93 | 0.73 | 0.60 | 0.28 |

As one can see, the efficiency is above 70% only up to 4 cores; above the
efficiency drops. Thus, for this system, state parallelization should only use
up to 4 cores.

Now the question is: why does the scaling break down above 4 processes and
especially above 8 processes? Let's look at the parallelization output from the
log of the 8-core run:

{{< code-block >}}
Info: Parallelization in states
Info: Node in group    0 will manage      4 states:     1 -      4
Info: Node in group    1 will manage      4 states:     5 -      8
Info: Node in group    2 will manage      4 states:     9 -     12
Info: Node in group    3 will manage      4 states:    13 -     16
Info: Node in group    4 will manage      4 states:    17 -     20
Info: Node in group    5 will manage      4 states:    21 -     24
Info: Node in group    6 will manage      4 states:    25 -     28
Info: Node in group    7 will manage      4 states:    29 -     32
{{< /code-block >}}

The system has 32 states, thus each core processes 4 states. In the previous
tutorial, it was indicated as a rule of thumb that 4 states per process is the
minimum to be efficient, also in terms of vectorization. As one can see, the
calculation time does not even decrease when going to 16 cores, where each
core processes only 2 states. However, on 4 cores, each process has 8 states and
that is more efficient.

##  Combined states and domain parallelization

Now run the system again, combining states and domain parallelization, with 2
cores each and with 2 and 4 cores. Extract the timings and compute the speed-up
and parallel efficiency. What is the most efficient way to run this system?

##  K-point parallelization

We don't treat k-point parallelization in detail here, because it is only
relevant for solids which are covered in a later tutorial. But as the
parallelization is quite trivial, you can scale the processes up to using only 1
k point per process. It is most efficient when the distribution of k points to
processes is balanced (e.g. all cores have 2 k-points instead of some having 2
and some having only 1).

##  Production runs

Usually, production runs are larger and require more resources. Thus, one
usually requests full nodes and then the scaling analysis is done using 1, 2, 4,
8, ... full nodes and the speed-up is computed relative to one node. Other than
that, the analysis is the same as outlined above.

For one of the systems you use in production runs: how many cores or nodes do
you usually use? Is that efficient? Run a quick scaling analysis on 3 or 4
different node counts to estimate the parallel efficiency. Check the guidelines
for the different parallelization strategies.


##  Memory usage

Another reason to use more nodes besides the larger compute power is the larger
memory that might be needed. If the memory needed for the states is too large to
fit into the main memory of one node, more nodes need to be used to distribute
the storage of the states across nodes.

There is a section in the output that gives approximate memory requirements. For
the example from the previous sections, it looks as follows:

{{< code-block >}}
#include_file doc/tutorials/hpc/scaling/1.domain_parallelization/memory.txt
{{< /code-block >}}

This indicates that mesh object takes about 10 MiB on each core and that the
states in total require 94 MiB for complex numbers (which are used for TD runs).
To estimate the total amount of memory needed, one can add the memory for the
states to the memory of mesh times the number of processes because the mesh
information is needed on each process. If this is larger than the memory per
node (which can usually be found in the documentation of the corresponding
cluster or supercomputer), one can estimate the number of nodes needed by
dividing the total memory needed by the memory per node.

Sometimes, supercomputers also offer compute nodes with more memory per node,
then it can also be an option to request those nodes.

##  Guidelines

Here is a brief summary of guidelines:
* do a scaling run for new systems and compute parallel efficiency
* parallel efficiency should be >70%
* in doubt use less cores and run several simulations in parallel
* for states parallelization: minimum 4-8 states per core
* for domain parallelization: ratio ghost/local points should be <25%
* for k-point parallelization: up to 1 k point per core
* for k-point and states parallelization: a balanced distribution is more efficient
* do not waste resources; leads to more science being done by all!

##  Code snippets

You can use the following python script to create a scaling plot as a log-log
plot of speed-up vs. core number:

```python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
# number of cores
cores = np.array([1, 2, 4, 8, 16])
# enter the times you measured here
times = np.array([0.90, 0.57, 0.34, 0.21, 0.11])
speedup = times[0]/times
ideal_speedup = cores/cores[0]
plt.loglog(cores, speedup, label="speed-up")
plt.loglog(cores, ideal_speedup, label="ideal")
plt.xlabel("Cores")
plt.ylabel("Speed-up")
plt.legend()
# some formatting
plt.gca().set_xticks(cores)
plt.gca().set_xticklabels(cores)
plt.gca().set_yticks(cores)
plt.gca().set_yticklabels(cores)
plt.gca().xaxis.set_tick_params(which='minor', size=0)
plt.gca().yaxis.set_tick_params(which='minor', size=0)
plt.savefig("speedup.pdf")
```

{{<tutorial-footer>}}
