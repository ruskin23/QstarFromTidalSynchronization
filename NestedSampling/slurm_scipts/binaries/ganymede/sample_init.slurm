#!/bin/bash

#SBATCH -J SampleInit
#SBATCH -n 15  #Total number of tasks to run at one time
#SBATCH -N 1  #Number of nodes on which to distribute the tasks evenly
#SBATCH -o /home/rxp163130/scratch/Nested_Sampling/slurm_output.76.1.out
#SBATCH -t 10:00:00

module load gsl

ulimit -S -m $((2*1024*1024))

export LD_LIBRARY_PATH=~/lib

MAIN_DIR="/home/rxp163130/QstarFromTidalSynchronization/NestedSampling/main"
OUT_DIR="/home/rxp163130/scratch/Nested_Sampling"

python3 $MAIN_DIR/main.py -s start -l 76 -n 15 >$OUT_DIR/sample_init.76.1.out
