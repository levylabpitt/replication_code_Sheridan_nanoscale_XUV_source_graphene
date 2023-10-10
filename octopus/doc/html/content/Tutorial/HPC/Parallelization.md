---
title: "Parallelization in octopus"
weight: 3
section: "Tutorial Parallelization_in_octopus"
series: "Tutorials"
tutorials: ["HPC"]
author: "Sebastian Ohlmann"
description: "Learn about the different strategies."
difficulties: "basic"
theories: "DFT"
calculation_modes: ["Ground state", "Time-dependent"]
system_types: "Molecule"
species_types: "Pseudopotentials"
features: "Total energy"


---


Octopus can run calculations in parallel to leverage the computing power of
clusters and it can utilize GPUs as well, but this is covered in another
tutorial. Octopus employs a hybrid parallelization of MPI and OpenMP to speed up
calculations. Using MPI implies that the code follows the ''distributed-memory''
paradigm: the data the code works on (in our case the Kohn-Sham wavefunctions)
are distributed over the processes and each process works on its local part of
the data. The OpenMP parallelization follows the ''shared-memory'' paradigm,
where work is distributed among threads that can access the same data.

In addition to the parallelization, Octopus employs vectorization which is a
capability of modern CPUs to execute several floating point operations
(additions and multiplications) in one instruction. How the data structures are
designed to make use of vectorization is part of an advanced tutorial.


##  Parallelization strategies

To distribute the data and the computations with MPI, Octopus follows different strategies that are orthogonal to each other and that can be combined:
* k points: for periodic systems, calculations usually comprise several k points to cover the Brillouin zone. Because the computation for each k point is mostly independent (except for computing the density), the parallelization scheme is simple and efficient: the k points are distributed over the processes and each process works on the local k points only.
* Kohn-Sham states: many operations in Octopus are the same for each Kohn-Sham state, especially for TD calculations. Thus this strategy is also quite efficient, but less so than the k point strategy because of the stronger coupling. Examples for the stronger coupling include orthogonalization and subspace diagonalization for GS runs.
* Domain: for this strategy, the real-space grid is distributed over the processes with a domain decomposition given by the METIS library. This strategy is more complicated and also less efficient than the other strategies, but needed when scaling to larger numbers of processes. The distribution of points leads to a certain number of inner points that each process works on. As Octopus uses finite difference stencils for computing derivatives, all points belonging to the stencil of inner points need to be available locally. When these points correspond to inner points on other processes, they are called ghost points and need to be updated before computing derivatives. They can also belong to the boundary, then their treatment depends on the boundary conditions. For this scheme to be efficient, the number of inner points must be large enough compared to the ghost points; otherwise the time for communicating the ghost points becomes too large.
* Electron-hole pairs: this strategy can be used for Casida calculations.

In addition to the MPI strategies, OpenMP can be used to exploit ''shared memory'' parallelization. 
When using OpenMP, each MPI process spawns a certain
number of OpenMP threads to parallelize computations and all threads have access
to the same memory of the process. In Octopus, this is used to parallelize loops
over the grid points, so it can be an alternative or addition to the
parallelization in domains.


##  First parallel runs

###  Ground state

For our first tests with parallel runs, we use the following input file, similar
to the time-dependent propagation tutorial:

{{< code-block >}}
#include_input doc/tutorials/hpc/parallelization/1.gs/inp
{{< /code-block >}}

Please use the batch script at the bottom of the tutorial to run octopus on one
core. You will get the output as described in the basic tutorials. Now modify
the batch script and submit it again to run on two cores. In case you run this
on your local PC, simply start octopus with `mpiexec -np <N> octopus`, where `<N>` is the number of cores to be used.

In the output, you will now see a new section on parallelization:

{{< code-block >}}
#include_input doc/tutorials/hpc/parallelization/1.gs/parallelization.txt
{{< /code-block >}}

As you can see, Octopus tells us that it was run with two processes and that it
used parallelization in domains ("ParDomains"). Below that, there is some
information on the parallelization in ions and then, more important, some more
information on the domain parallelization:

{{< code-block >}}
#include_input doc/tutorials/hpc/parallelization/1.gs/parmetis.txt
{{< /code-block >}}

First, the output tells you how the partitioning of the mesh was done, in this
case using ParMETIS, the parallel version of METIS. Very importantly, you can
see the information on the mesh partitioning: for each process ("node"), it
shows the number of neighbours, and the number of local, ghost, and boundary
points.

This information is useful to judge the efficiency of the domain
parallelization: as a rule of thumb, if the number of ghost points approaches
25% of the number of local points, the scheme usually becomes inefficient.

The number of boundary points only depends on the size of the stencil and can
usually not be changed. For finite systems, their impact is small; for periodic
systems they need to be updated which can impact the performance.


###  Time-dependent calculation

Now run a time-dependent calculation for the same system using the following
input file:

{{< code-block >}}
#include_input doc/tutorials/hpc/parallelization/2.td/inp
{{< /code-block >}}


You can use the same batch script as for the ground-state run to execute octopus
on two cores.

When you look at the output, you should see the parallelization section as
follows:

{{< code-block >}}
#include_input doc/tutorials/hpc/parallelization/2.td/parallelization.txt
{{< /code-block >}}

As you can see, octopus has used 2 processes as expected. Both processes are
used for parallelization in states ("ParStates"). Slightly further down you see
how the states are distributed:

{{< code-block >}}
#include_input doc/tutorials/hpc/parallelization/2.td/parstates.txt
{{< /code-block >}}

This means that the both processes have 2 states that they manage.

Because this system is quite small, you will probably not see a big change in
run time for the different parallelization options.


##  Control parallelization strategies

By default, octopus chooses a certain combination of parallelization
strategies. However, this can (and should) be changed by the user who knows best
about the system and also about the machine where the code is executed.

The following variables control the parallelization strategies:
* {{<variable "ParKPoints">}}: parallelization in k-points
* {{<variable "ParStates">}}: parallelization in Kohn-Sham states
* {{<variable "ParDomains">}}: parallelization in domains, i.e. for grid points
* {{<variable "ParOther">}}: "other" parallelization; mostly used for electron-hole pairs in the Casida mode
These variables can be set to the number of processes used for each of the
strategies or to "auto" to let octopus decide about the number of processes or
to "no" to disable the corresponding strategy. By default, all variables are set
to "auto" for TD calculations. For all other modes, ParStates is set to "no",
whereas the others are still set to "auto".

When you set the number of processes for the different parallelization
strategies, '''their product must be equal to the total number of processes'''
used in this computation. This is due to the fact that the strategies are
orthogonal to each other. Using different strategies corresponds to dividing a
hypercube in different dimensions.

Let's run the TD example from above again using two processes, but this time
with domain parallelization. For this, you need to add `ParDomains = 2` to the
input file and submit the job script to run on two processes. The
parallelization section in the output will then look like

{{< code-block >}}
#include_input doc/tutorials/hpc/parallelization/3.strategies/parallelization.txt
{{< /code-block >}}

This shows you that indeed the parallelization in domains has been used. A bit
further down you can again check how the mesh is distributed (which should be
the same as in the GS run).

Now you can combine parallelization in states and domains: add `ParStates = 2`
to the input file and change the slurm batch script to start 4 processes (i.e.
use `--ntasks=4`). Then submit the job. The parallelization section in the
output will be:

{{< code-block >}}
#include_input doc/tutorials/hpc/parallelization/3.strategies/parstates.txt
{{< /code-block >}}

From this output you can see that octopus used 2 processes for states
parallelization and 2 processes for domain parallelization. This means that each
process will handle 2 states and about half of the grid.


##  OpenMP parallelization

To try out the OpenMP parallelization, you can use the same input file as
before, but delete all the parallelization variables. To use OpenMP, get the
correct batch script (see at the end of the tutorial), which sets
`OMP_NUM_THREADS` environment variable. If you run octopus on your laptop, you need
to set this variable to the number of OpenMP threads you would like to start.

If you run this, you will see the following information on the parallelization:

{{< code-block >}}
#include_input doc/tutorials/hpc/parallelization/4.omp/parallelization.txt
{{< /code-block >}}

This indicates that one MPI process is used with two OpenMP threads per process.

The OpenMP parallelization can be combined with the MPI parallelization - their
product must be equal to the total number of CPU cores used by the batch job.

Now you can repeat the last example from the previous section, but in addition
use 2 OpenMP threads (so add again the Par* variables and set `--ntasks=4`). The
output will show:

{{< code-block >}}
#include_input doc/tutorials/hpc/parallelization/5.hybrid/parallelization.txt
{{< /code-block >}}

As expected, there are 4 processes, 2 for states and 2 for domains
parallelization and there are 2 OpenMP threads per MPI process. In total, this
run uses 8 CPU cores.


##  Compare performance of different strategies

Let us compare the performance of the different strategies for a different
molecule that is slightly larger. Save the following as inp file:

{{< code-block >}}
#include_input doc/tutorials/hpc/parallelization/6.compare/inp
{{< /code-block >}}

and the following as "1ala.xyz":

{{< code-block >}}
#include_file doc/tutorials/hpc/parallelization/6.compare/1ala.xyz
{{< /code-block >}}

Then, run the ground state calculation.

Now, run the TD calculation (change the CalculationMode variable to td) using one
process. After the initialization you see something like the following output:

{{< code-block >}}
  *  *  *  *  *  *  *  *  *  ** Time-Dependent Simulation   *  *  *  *  *  *  *  *  *  *  *
  Iter           Time        Energy   SC Steps    Elapsed Time

  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *    *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
      1       0.050000   -109.463446         1         1.087
      2       0.100000   -109.463446         1         1.076
      3       0.150000   -109.463446         1         1.077
      4       0.200000   -109.463446         1         1.089
      5       0.250000   -109.463446         1         1.082
      6       0.300000   -109.463446         1         1.064
      7       0.350000   -109.463446         1         1.070
      8       0.400000   -109.463446         1         1.100
      9       0.450000   -109.463446         1         1.086
     10       0.500000   -109.463446         1         1.186
[...]
{{< /code-block >}}

The last column indicates how long one timestep took in real time, in this case
about 1.1 s.

To assess the performance of the different parallelization strategies, use your
knowledge from the previous sections and run the TD calculation:
* on two cores using states parallelization
* on two cores using domain parallelization
* on one core and two OpenMP threads
For each of these, get an average time per timestep from the output.

Compare the time of each of those runs with the serial one. How big is the
speed-up? The speed-up is defined as the time for the serial run divided by the
time for the run on two processes. Ideally this number should be two. How far
away from that are your runs? Which mode is most efficient in this case?

For the domain-parallel run: what is the ratio of ghost points to local points?
Is it still ok?

For the states-parallel run: how many states per process do you have? Is that
still ok? Is the distribution balanced?

If you still have time, compare the following runs on 4 cores in the same way:
* 4 cores using states parallelization
* 4 cores using domain parallelization
* 2 cores using states parallelization + 2 cores using domain parallelization
and answer the same questions to check if the runs are efficient.


##  Tips and tricks

If you run on a cluster, make sure to read the documentation to understand how
many cores per node it has. Try to use the full number of processes available on
the number of nodes you request. Do not use hyperthreading, octopus will not
profit from this.

For k-point-parallelization, you can go down to one k point per process. For
states parallelization, you should still have 4 states per process as a minimum
to effectively use vectorization. For both of these modes, the best distribution
is balanced, i.e., each process has the same number of k points or states.

For domain parallelization, usually the number of ghost points should be at
maximum 25% of the local points. For large grids, it can be more efficient to
use up to 10 OpenMP threads instead of processes in the domain parallelization.
As an example, if you use `ParDomains = 32`, you can instead try `ParDomains =
8` and set `OMP_NUM_THREADS=4` to run with 4 OpenMP threads per process. Be
aware that the number of OpenMP threads should be smaller than the number of
cores per socket in a node to avoid problems with non-uniform memory access
(NUMA). For small grids (< about 10^5 points), it is usually not worth to use
either OpenMP or domain parallelization.

It is often worth trying to find a good factorization of the total number of
processes to use the different parallelization strategies efficiently.


##  Slurm batch scripts for MPCDF systems

To run on one core of the reserved resources, please use the following batch
script:

```bash
#include_file doc/tutorials/hpc/parallelization/7.scripts/job_mpi.sh
```

To run on two cores, please replace `--ntasks=1` by `--ntasks=2`.

To utilize OpenMP, please use the following batch script:
```bash
#include_file doc/tutorials/hpc/parallelization/7.scripts/job_omp.sh
```

This will run with one MPI rank and 2 OpenMP threads. To change the number of
OpenMP threads per MPI rank, adapt `--cpus-per-task` accordingly.

{{<tutorial-footer>}}