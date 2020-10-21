#!/bin/bash

#SBATCH -J MCMCU
#SBATCH -p debug 
#SBATCH -n 1  #Total number of tasks to run at one time
#SBATCH -N 1  #Number of nodes on which to distribute the tasks evenly
#SBATCH -o /home/rxp163130/scratch/debug/debug.out
#SBATCH -t 00:20:00

module load launcher 
module load gsl 

ulimit -S -m $((2*1024*1024))

export LD_LIBRARY_PATH=~/lib

python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/version1/main.py -c -i 5 -l 93 -m uncorrelated >/home/rxp163130/scratch/Run28/93.5.out

