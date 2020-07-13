#!/bin/bash

#SBATCH -J MCMCW 
#SBATCH -n 185  #Total number of tasks to run at one time
#SBATCH -N 12
#SBATCH -o home/rxp163130/scratch/version2/first_run/slurm.output.out
#SBATCH -t 4-00:00:00

module load launcher 
module load gsl

ulimit -S -m $((2*1024*1024))

export LD_LIBRARY_PATH=~/lib
export LAUNCHER_JOB_FILE=jobfile
 
$LAUNCHER_DIR/paramrun


