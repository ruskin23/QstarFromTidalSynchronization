#!/bin/bash

#SBATCH -J binaries 
#SBATCH -n 81  #Total number of tasks to run at one time
#SBATCH -N 6  #Number of nodes on which to distribute the tasks evenly
#SBATCH -o /home/rxp163130/scratch/bianry_evolution_analysis/binaries.out
#SBATCH -t 2-00:00:00

module load launcher 
module load gsl

ulimit -S -m $((2*1024*1024))

export LD_LIBRARY_PATH=~/lib
export LAUNCHER_JOB_FILE=jobfile
 
$LAUNCHER_DIR/paramrun


