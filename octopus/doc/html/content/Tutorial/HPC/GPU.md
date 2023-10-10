---
title: "Octopus on GPUs"
weight: 5
section: "Tutorial Octopus_on_GPUs"
series: "Tutorials"
tutorials: ["HPC"]
author: "Sebastian Ohlmann"
description: "Utilizing the GPUs."
---


Octopus supports using GPUs and this support has been substantially improved in
the last years. It does not yet support all features, but for certain systems,
using GPUs can be very efficient. In this tutorial, you will learn how to run
Octopus on GPUs and also for which systems you can expect the best efficiency.

##  Supported Features

In general, the most-often used features should work on GPUs and should be efficient as
well. If you run ground-state calculations, please use the RMMDIIS eigensolver. Time-dependent
calculations are supported and are the most efficient (and also the most-used) run mode.
Some features do not work when running on GPUs or are quite inefficient, e.g. spin-orbit coupling, DFT+U, hybrid functionals.
Moreover, there can be an issue for periodic systems when spheres of pseudopotentials overlap with themselves due to the periodic boundary conditions.

It is really worth trying your production workloads on GPUs - for quite a lot of applications, it is more
efficient than on CPUs and there will be more and more GPU resources available in the future.

If you find any issues with functionality or performance, please notify the developers.

##  Implementation

The implementation is based on CUDA and targets only NVIDIA GPUs at the moment.
As Octopus needs double-precision floating point operations, most consumer-grade
and gaming GPUs will only deliver inferior performance because they focus on
single-precision operations. The high-end HPC GPUs that are usually available in
computing centers deliver much better double-precision performance.

When octopus is compiled with GPU support, it will detect if a GPU is available.
If there are several GPUs per node, you should generally use one MPI process per
GPU to be most efficient.

In general, systems need to be large enough in number of grid points and number
of states to expose enough parallelism that the GPUs can exploit. If the system
is too small, the overhead of launching computations on the GPUs and of copying
data from and to the GPUs makes the computations less efficient than on CPUs.
For large systems, the performance of running on GPUs can be much larger than
the electricity consumption and the costs of the GPUs compared to CPUs.

##  Run a calculation on the GPU

For this tutorial, we will use the same input as in the {{<tutorial "HPC/Scaling" "previous tutorial">}}.
Save the following as {{< file "inp" >}} file:

{{< code-block >}}
#include_input doc/tutorials/hpc/gpu/inp
{{< /code-block >}}

and the following as {{< file "1ala.xyz" >}}:

{{< code-block >}}
#include_input doc/tutorials/hpc/gpu/1ala.xyz
{{< /code-block >}}

Then, run the ground state calculation or reuse the ground state from the
previous tutorial.

You can now run your first calculation using a GPU: change the run mode to td
and use the batch script to run on one GPU as given at the end of this page.
In the output, you will see the following section that shows some information on
the GPU being used:

{{< code-block >}}
************************** GPU acceleration **************************
Rank 00000 uses device number 00000

Selected device:
      Framework              : CUDA
      Device type            : GPU
      Device vendor          : NVIDIA Corporation
      Device name            : Tesla V100-PCIE-32GB
      Cuda capabilities      : 7.0
      Driver version         : 11000
      Device memory          : 32510.500 MiB
      Local/shared memory    : 48.000 KiB
      Max. group/block size  : 1024

**********************************************************************
{{< /code-block >}}

The parallelization section shows that one process is used because we want one
process per GPU. Because we use half a node, we can still use 20 OpenMP threads
which can be beneficial especially in the initialization where still some
computations are done on CPUs.

At the end of the run, there is an additional section on allocations and
deallocations that shows how well the internal caching is used. Normally, you
can ignore this section.

Now extract the timing of one timestep as we did in the previous tutorial. As a
reference, I have obtained 0.074s. To compare this using the same resources, but
no GPU, you can set the input variable {{<variable "DisableAccel">}} {{<code " = yes">}} and submit
the same batch script again. This time, it will only run on CPUs. For this, I
get a time of 0.26. We can compute a speed-up from these timings, which is 3.5;
for your runs it might deviate slightly.

This means that we get a speed-up factor of 3.5 when comparing half a GPU node
on cobra with half a CPU node for this rather small system. This is already
nice!

##  Beneficial use cases

TD calculations are most optimized for GPUs, other run modes less so. For GS
runs, e.g., the orthogonalization and subspace diagonalization cannot make full
use of GPUs yet, which leads to smaller benefits for GS runs on GPUs.

Moreover, grids need to be large enough to saturate the GPUs which are very
powerful. Each process should have more than $10^5$ grid points.

NVIDIA GPUs are built such that 32 threads always execute the same instructions.
Because of this, it is most efficient to run octopus in such a way that each
process has a multiple of 32 states. The internal data structures of the code
are designed to expose this parallelism in the threads to the GPU.

Parallelization in states (and k points) is more efficient than parallelization
in domains. This is also the case for CPU runs, but the different is more
pronounced on GPUs because the communication of the boundary points necessitates
also transfers between the GPUs or between the GPU and the CPU. We use
CUDA-aware MPI to enable direct data transfers between the GPUs, but this cannot
completely compensate the overhead.

##  Multi-GPU runs

To examine how octopus runs using several GPUs, you can use the second jobscript
at the end of the tutorial to submit a job to a full GPU node, using 2 MPI
processes and 2 GPUs. The output section on GPUs will show:
{{< code-block >}}
************************** GPU acceleration **************************
Rank 00001 uses device number 00001
Rank 00000 uses device number 00000

Selected device:
      Framework              : CUDA
      Device type            : GPU
      Device vendor          : NVIDIA Corporation
      Device name            : Tesla V100-PCIE-32GB
      Cuda capabilities      : 7.0
      Driver version         : 11000
      Device memory          : 32510.500 MiB
      Local/shared memory    : 48.000 KiB
      Max. group/block size  : 1024

**********************************************************************
{{< /code-block >}}

which confirms that 2 GPUs are used. As you can see in the parallelization
section, by default states parallelization is used, which results in 16 states
per process.

Extract the timing from the output. I obtain 0.043s, corresponding to a speed-up
of about 1.72 and a parallel efficiency of 86%. Although 16 states per process
is not optimal, the parallel efficiency is still very good.

Now run the same calculation on 2 full nodes and determine the timing and
compute the speed-up and parallel efficiency. What do you see? Is it still
efficient?

Next, we will run the same calculation on one node, but using domain
parallelization. For this, add {{<variable "ParStates">}}<tt> = 1</tt> to the input file and
submit the jobscript again. Examine the parallelization section in the output to
ensure that domain parallelization has been used. Extract the timing from the
output and compute the speed-up and parallel efficiency. I obtain about 0.090s
which gives a speed-up of 0.82 - the code becomes slower! This is due to the
small grid in this system which does not expose enough parallelism to distribute
it efficiently to two GPUs - the overhead of the parallelization scheme (e.g.
communicating ghost points) is too large.

##  Domain parallelization for larger grids

Set the spacing to `0.15` and rerun the ground state calculation (you
need to delete the restart folder before). Then repeat the runs on one GPU and
on 2 GPUs using parallelization in domains. What do you observe? Is this
parallelization scheme more efficient for the larger grid?


##  Compiling Octopus for GPUs

If you compile octopus yourself, you need to specify the flag
`--enable-cuda` to enable CUDA support. If your CUDA installation is in a
non-standard path, you additionally need to specify `--with-cuda-prefix=DIR` to
point to the corresponding CUDA installation directory. If you would like to do
profiling on GPUs using Nsight systems, it is advisable to add the flag `--enable-nvtx` which will add
support for NVTX.


##  Using NVLink on Raven

If the domain parallelization on larger grids is not efficient on cobra, you can try to use a system that features
a fast interconnect (e.g. NVLink), such as the raven supercomputer at MPCDF. If you have a normal account at MPCDF: run the same
example on raven. The GPU nodes on raven have a fast interconnect on the board
(NVLink) which can massively speed up domain parallelization runs which use
CUDA-aware MPI. You need to load the modules "octopus-gpu/12" in order to use the 
octopus version specifically built to make
use of CUDA-aware MPI. Run the example on 1, 2, and 4 GPUs with domain
parallelization. What do you observe? Is it more efficient on raven?



##  Slurm scripts

To run octopus on the cobra supercomputer using one GPU, you can use the
following script:

```bash
#!/bin/bash -l
# Standard output and error:
#SBATCH -o ./tjob.out.%j
#SBATCH -e ./tjob.err.%j
# Initial working directory:
#SBATCH -D ./
#
# Job Name:
#SBATCH -J octopus_course
#
# Reservation:
#SBATCH --reservation=mpsd_course
#
# Node feature:
#SBATCH --constraint="gpu"
# Specify type and number of GPUs to use:
#SBATCH --gres=gpu:v100:1       - Use only 1 GPU of a shared node
#SBATCH --mem=92500             - Memory is necessary if using only 1 GPU
#
# Number of nodes and MPI tasks per node:
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
# for OpenMP:
#SBATCH --cpus-per-task=20
#
#SBATCH --mail-type=none
#SBATCH --mail-user=userid@example.mpg.de
#
# wall clock limit:
#SBATCH --time=00:05:00

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
# For pinning threads correctly:
export OMP_PLACES=cores

# Run the program:
module purge
module load octopus-gpu/{{<octopus-version>}}
srun octopus
```

It will run on a node with NVIDIA V100 GPUs and use 1 of them.

To run on a full cobra node with 2 GPUs, you can use the following script:

```bash
#!/bin/bash -l
# Standard output and error:
#SBATCH -o ./tjob.out.%j
#SBATCH -e ./tjob.err.%j
# Initial working directory:
#SBATCH -D ./
#
# Job Name:
#SBATCH -J octopus_course
#
# Reservation:
#SBATCH --reservation=mpsd_course
#
# Node feature:
#SBATCH --constraint="gpu"
# Specify type and number of GPUs to use:
#SBATCH --gres=gpu:v100:2       - Use both GPUs of a node
#
# Number of nodes and MPI tasks per node:
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
# for OpenMP:
#SBATCH --cpus-per-task=20
#
#SBATCH --mail-type=none
#SBATCH --mail-user=userid@example.mpg.de
#
# wall clock limit:
#SBATCH --time=00:05:00

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
# For pinning threads correctly:
export OMP_PLACES=cores

# Run the program:
module purge
module load octopus-gpu/{{< octopus-version >}}
srun octopus
```

{{< tutorial-footer >}}
