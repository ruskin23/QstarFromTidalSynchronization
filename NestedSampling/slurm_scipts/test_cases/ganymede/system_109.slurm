#!/bin/bash

#SBATCH -J S109n4
#SBATCH -n 4
#SBATCH -N 1
#SBATCH -t 4-00:00:00
#SBATCH -o /home/rxp163130/scratch/Nested_Sampling/System109_slurm_n4.out

module load gsl

ulimit -S -m $((2*1024*1024))

export LD_LIBRARY_PATH=~/lib

python3 /home/rxp163130/QstarFromTidalSynchronization/NestedSampling/slurm_scipts/main.py -l 109 >/home/rxp163130/scratch/Nested_Sampling/109.n4.out
