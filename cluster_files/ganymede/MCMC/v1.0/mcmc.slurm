#!/bin/bash

#SBATCH -J MCMC 
#SBATCH -n 155 #Total number of tasks to run at one time
#SBATCH -N 10  #Number of nodes on which to distribute the tasks evenly
#SBATCH -o /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/output/ganymede/MCMC.o%j
#SBATCH -t 4-00:00:00

module load launcher 
module load gsl

export LD_LIBRARY_PATH=~/lib
export LAUNCHER_JOB_FILE=jobfile
 
$LAUNCHER_DIR/paramrun


