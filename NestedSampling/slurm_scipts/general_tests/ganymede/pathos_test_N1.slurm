#!/bin/bash

#SBATCH -J NestedN1
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 06:30:00
#SBATCH -o /home/rxp163130/scratch/Nested_Sampling/slurm_output_N1.out

module load gsl

ulimit -S -m $((2*1024*1024))

export LD_LIBRARY_PATH=~/lib

python3 /home/rxp163130/QstarFromTidalSynchronization/NestedSampling/slurm_scipts/pathos_test_main.py -l 1 -n 1 >/home/rxp163130/scratch/Nested_Sampling/1.N1.out
