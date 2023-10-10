#!/bin/bash -l
# Standard output and error:
#SBATCH -o ./tjob.out.%j
#SBATCH -e ./tjob.err.%j
# Initial working directory:
#SBATCH -D ./
# Job Name:
#SBATCH -J td-gr-TWcm2-501
# Queue (Partition):
#SBATCH --partition=mpsd
# Number of nodes and MPI tasks per node:
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#
#SBATCH --mail-type=ALL
#SBATCH --mail-user=davis.welakuh@mpsd.mpg.de
#
# Wall clock limit:
#SBATCH --time=12:00:00
#
# For pinning threads correctly:
export OMP_PLACES=cores
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

module purge
module load intel/19.1.3 impi/2019.9 mkl/2020.1 gsl libxc/5.1.5 hdf5-serial netcdf-serial cuda metis parmetis fftw-mpi
module list

# Run the program:
srun $HOME/td_photon/latest_octopus/local_build/installed/bin/octopus > td.out
