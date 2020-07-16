#!/bin/bash

#SBATCH -J S109n12
#SBATCH -n 12
#SBATCH -N 1
#SBATCH -t 03:30:00
#SBATCH -o /home/rxp163130/scratch/Nested_Sampling/System109_slurm_n12.out

module load gsl

ulimit -S -m $((2*1024*1024))

export LD_LIBRARY_PATH=~/lib

python3 /home/rxp163130/QstarFromTidalSynchronization/NestedSampling/slurm_scipts/main.py -l 109 -n 12 >/home/rxp163130/scratch/Nested_Sampling/109.n12.out
