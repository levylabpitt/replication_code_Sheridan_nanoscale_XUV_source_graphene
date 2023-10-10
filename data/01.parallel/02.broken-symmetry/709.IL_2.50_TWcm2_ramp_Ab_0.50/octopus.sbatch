#!/bin/bash -l
# Standard output and error:
#SBATCH -o ./tjob.out.%j
#SBATCH -e ./tjob.err.%j
# Initial working directory:
#SBATCH -D ./
# Job Name:
#SBATCH -J 006-td-gr-TWcm2
# Queue (Partition):
##SBATCH --partition=mpsd
# Number of nodes and MPI tasks per node:
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=40
#SBATCH --cpus-per-task=1
#
#SBATCH --mail-type=ALL
#SBATCH --mail-user=davis.welakuh@mpsd.mpg.de
#
# Wall clock limit:
#SBATCH --time=10:00:00
#
# For pinning threads correctly:
export OMP_PLACES=cores
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

module purge
module load gcc/11 openmpi/4
module list

# Run the program:
srun $HOME/TDDFT-real-time/octopus/local_build/installed/bin/octopus > td.out
