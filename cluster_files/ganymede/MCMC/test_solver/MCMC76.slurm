#!/bin/bash

#SBATCH -J MCMC76 
#SBATCH -n 10  #Total number of tasks to run at one time
#SBATCH -N 1  #Number of nodes on which to distribute the tasks evenly
#SBATCH -o /home/rxp163130/scratch/test_solver/slurm_output_76.out
#SBATCH -t 4-00:00:00

module load launcher 
module load gsl

ulimit -S -m $((2*1024*1024))

export LD_LIBRARY_PATH=~/lib
export LAUNCHER_JOB_FILE=jobfile
 
$LAUNCHER_DIR/paramrun


