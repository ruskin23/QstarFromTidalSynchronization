#!/bin/bash

#SBATCH -J NestedN4
#SBATCH -n 4
#SBATCH -N 4
#SBATCH -t 00:30:00
#SBATCH -o /home/rxp163130/scratch/Nested_Sampling/slurm_output.out

module load gsl

ulimit -S -m $((2*1024*1024))

export LD_LIBRARY_PATH=~/lib

python3 /home/rxp163130/QstarFromTidalSynchronization/NestedSampling/slurm_scipts/pathos_test_main.py -l 1 -n 4 >/home/rxp163130/scratch/Nested_Sampling/1.n4.out
