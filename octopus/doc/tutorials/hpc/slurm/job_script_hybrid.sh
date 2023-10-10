#!/bin/bash -l
# Standard output and error:
#SBATCH -o ./tjob.out.%j
#SBATCH -e ./tjob.err.%j
# Initial working directory:
#SBATCH -D ./
# Job Name:
#SBATCH -J octopus_course
#
# Reservation:
#SBATCH --reservation=mpsd_course
#
#
# Number of nodes and MPI tasks per node:
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
# for OpenMP:
#SBATCH --cpus-per-task=4
#
#SBATCH --mail-type=none
#SBATCH --mail-user=userid@example.mpg.de
#
# Wall clock limit:
#SBATCH --time=00:01:00

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
# For pinning threads correctly:
export OMP_PLACES=cores

# Run the program:
module purge
module load octopus/{{<octopus-version>}}
srun octopus

