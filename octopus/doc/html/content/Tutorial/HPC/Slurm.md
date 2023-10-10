---
title: "Slurm usage"
weight: 2
section: "Tutorial Slurm_usage"
series: "Tutorials"
tutorials: ["HPC"]
author: "Sebastian Ohlmann"
description: "The job submission system."
---


##  HPC systems

HPC systems provide large complex shared resources to many users in parallel.
Therefore, central workload management is mandatory for users and admins.
This is typically achieved by using a resource manager and a job scheduler.
The resource manager knows about all the resources available on a system
and monitors the availability and load of each node. It basically manages all
resources like CPUs, memory, and GPUs in a cluster. The job scheduler
assigns compute tasks to resources and manages queues of compute tasks and
their priority. It is typically configured for best possible utilization of the
resources by given typical workloads.


##  Slurm

Slurm is a resource manager and job scheduler that is used by the majority of
the TOP 500 HPC systems and that is used on all HPC systems and clusters at
MPCDF. It is open source software with commercial support (Documentation:
https://slurm.schedmd.com, MPCDF HPC documentation:
https://docs.mpcdf.mpg.de/doc/computing/).

Some slurm terminology:
* Job: Reservation of resources on the system to run job steps
* Job step: program/command to be run within a job, initiated via `srun`
* Node: physical multi-core shared-memory computer, a cluster is composed of many nodes
* CPU: single processing unit (core), a node contains multiple CPUs
* Task: process (i.e. instance of a program being executed), may use one or more CPU up to all CPUs available a node, a job step may run multiple tasks in parallel over several nodes
* Partition: a “queue”, where to run jobs, defines specific resource limits or access control


##  Slurm commands

Important slurm commands are:
* `sinfo`: show state of partitions and resources managed by slurm
* `sbatch job_script.sh`: submit a job script for later execution, obtain job_id
* `scancel job_id`: cancel a job, or send signals to tasks
* `squeue`: show state of jobs or job steps in priority order
* `srun executable`: Initiate job step, launch executable (typically used in job scripts)
* `sacct`: show information for finished jobs

You can get a list of waiting and running jobs of yourself with `squeue --me`.

You can display a concise list of partitions with `sinfo -s` (A/I/O/T means
allocated/idle/offline/total).


##  Slurm jobs

Slurm jobs are submitted by users from the login node and then scheduled by
slurm to be executed on the compute nodes of the cluster.

Any slurm job requires:
- Specification of the resources – “what does the job need?”
  * Duration
  * Number of CPUs
  * Amount of memory
  * GPUs
  * other resources or constraints
- Definition of the job steps – “what should the job do?”
  * commands/programs to be executed via `srun`
  * typically, first some module are loaded, then the program is executed using srun

All this information is bundled in job scripts that are submitted to slurm
utilizing the `sbatch` command.


##  Submitting a first job script

Let's create a job script to run a simple octopus calculation.

First, generate the input file called {{< file "inp" >}} with a text editor (the same as
in the {{< tutorial "Basics/Getting_started" "very first tutorial">}}):

{{< code-block >}}
#include_input doc/tutorials/octopus_basics/getting_started/1.H_atom/inp
{{< /code-block >}}

Second, create a job script file called {{< file "job_script.sh" >}}:

```bash
#include_input doc/tutorials/hpc/slurm/job_script.sh
```

For this job, output will be written to {{< file "tjob.out.XXX" >}} (XXX is the job id), error
output to {{< file "tjob.err.XXX" >}}. The job will run using one MPI task on one core,
requesting 2200 MB of memory. You can change the `--mail-type` option to `all`
and the `--mail-user` option to your email address to get email notifications
about changes in the job status (i.e. when the job starts and ends). The job
requests a time of one minute (`--time` option). In the script, a octopus
module is loaded and the octopus executable is started with srun. The
`--reservation` option is only needed for this course to use dedicated
resources reserved for us.

Now submit the job with `sbatch job_script.sh`. You should see an output like

```
Submitted batch job XXX
```

where XXX is the job id. You can check the status by running `squeue --me`.

Once the calculation has finished, you can check the output by opening the file
{{< file "tjob.out.XXX" >}}. Moreover, you should see the folders and files that octopus has
created, as in the first tutorial.


##  More job script examples

To submit a parallel job on a few cores, but still on one node (cobra has 40
cores per node), you can use the options `--ntasks=8` and `--mem=17600` to run
on 8 cores, for example.

To run octopus on a full node (or several full nodes) in pure MPI mode, please
use the following job script:

```bash
#include_input doc/tutorials/hpc/slurm/job_script_mpi.sh
```

Save this file as {{<file "job_script_mpi.sh">}} and submit it with `sbatch
job_script_mpi.sh`.  This will run octopus on all 40 cores of one node. To run
on multiple nodes, adapt the `--nodes` option accordingly.

To run octopus in hybrid mode (MPI + OpenMP), which is suitable for large
grids, you can employ the following script:

```bash
#include_input doc/tutorials/hpc/slurm/job_script_hybrid.sh
```

This will run octopus on one full node, using 10 MPI ranks with 4 OpenMP
threads each. Exporting the environment variables is necessary to ensure
correct pinning of all processes and threads.

{{< tutorial-footer >}}
