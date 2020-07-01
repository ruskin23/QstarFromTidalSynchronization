#!/bin/bash

#SBATCH -J NestedTest
#SBATCH -n 4
#SBATCH -N 1
#SBATCH -t 1:00:00
#SBATCH -o /home/rxp163130/scratch/Nested_Sampling/slurm_output.out

module load gsl

ulimit -S -m $((2*1024*1024))

export LD_LIBRARY_PATH=~/lib

python3 /home/rxp163130/QstarFromTidalSynchronization/NestedSampling/slurm_scipts/main.py -l 111 >/home/rxp163130/scratch/Nested_Sampling/111.out
