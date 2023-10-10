#!/bin/bash -l
#SBATCH -J pulpo
#SBATCH -N 2
#SBATCH -p regular
#SBATCH -t 04:00:00
#SBATCH --export=ALL

cd $HOME/edison/octopus
export OMP_NUM_THREADS=1
export MPIEXEC=`which srun`
export TEMPDIRPATH=$SCRATCH/tmp
export OCT_TEST_NJOBS=15
make check &> $SLURM_SUBMIT_DIR/makecheck