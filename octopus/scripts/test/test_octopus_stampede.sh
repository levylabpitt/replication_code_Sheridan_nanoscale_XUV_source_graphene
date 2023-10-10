#!/usr/bin/env bash

#SBATCH -n 16
#SBATCH -p development
#SBATCH -t 01:00:00
#SBATCH --export=ALL

module load perl
WORKDIR=$PWD
export TEMPDIRPATH=$SCRATCH/tmp
cd $HOME/octopus
export OCT_TEST_NJOBS=8
export MPIEXEC=`which ibrun`
make check &> $WORKDIR/makecheck