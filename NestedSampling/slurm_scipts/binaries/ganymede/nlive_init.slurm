#!/bin/bash

#SBATCH -n 80  #Total number of tasks to run at one time
#SBATCH -N 5  #Number of nodes on which to distribute the tasks evenly
#SBATCH -t 1-00:00:00

module load gsl

ulimit -S -m $((2*1024*1024))

export LD_LIBRARY_PATH=~/lib

MAIN_DIR="/home/rxp163130/QstarFromTidalSynchronization/NestedSampling/main"
OUT_DIR="/home/rxp163130/scratch/Nested_Sampling"

python3 $MAIN_DIR/main.py -s start -l $S -n 80 >$OUT_DIR/nlive_init.$S.out
