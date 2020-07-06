#!/bin/bash

#SBATCH -J NestedN8
#SBATCH -n 8
#SBATCH -N 8
#SBATCH -t 00:30:00
#SBATCH -o /home/rxp163130/scratch/Nested_Sampling/slurm_output_N8.out

module load gsl

ulimit -S -m $((2*1024*1024))

export LD_LIBRARY_PATH=~/lib

python3 /home/rxp163130/QstarFromTidalSynchronization/NestedSampling/slurm_scipts/pathos_test_main.py -l 1 -n 8 >/home/rxp163130/scratch/Nested_Sampling/1.N8.out
