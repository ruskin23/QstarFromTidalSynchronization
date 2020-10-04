#!/bin/bash

#SBATCH -J MCMCU 
#SBATCH -n 65  #Total number of tasks to run at one time
#SBATCH -N 5  #Number of nodes on which to distribute the tasks evenly
#SBATCH -o /home/rxp163130/scratch/Run25/MCMCNonSync.o%j
#SBATCH -t 4-00:00:00

module load launcher 
module load gsl

ulimit -S -m $((2*1024*1024))

export LD_LIBRARY_PATH=~/lib
export LAUNCHER_JOB_FILE=uncorrelated
 
$LAUNCHER_DIR/paramrun


