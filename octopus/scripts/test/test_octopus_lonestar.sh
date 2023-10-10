#!/bin/bash

#$ -N test_pulpo
#$ -cwd
#$ -pe 6way 12
#$ -q development
#$ -l h_rt=00:45:00
#$ -V

cd $HOME/octopus
export TEMPDIRPATH=$SCRATCH/tmp
if [ ! -d $TEMPDIRPATH ]; then
    mkdir $TEMPDIRPATH
fi

export MPIEXEC=`which ibrun`
export OCT_TEST_NPROCS=1
make check &> makecheck