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
# Number of MPI Tasks, e.g. 1:
#SBATCH --ntasks=1
#SBATCH --ntasks-per-core=1
# Memory usage [MB] of the job is required, 2200 MB per task:
#SBATCH --mem=2200
#
#SBATCH --mail-type=none
#SBATCH --mail-user=userid@example.mpg.de
#
# Wall clock limit:
#SBATCH --time=00:05:00

# Run the program:
module purge
module load octopus/{{<octopus-version>}}
srun octopus

