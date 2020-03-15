#!/bin/bash

#SBATCH -J PercentileAges2 
#SBATCH -n 128 #Total number of tasks to run at one time
#SBATCH -N 8  #Number of nodes on which to distribute the tasks evenly
#SBATCH -o /home/rxp163130/QstarFromTidalSynchronization/binary_star_evolution/analyze_spin_v_logQ/general_spin_v_logQ/output/ganymede/PercentileAges/PercentileAges.o%j
#SBATCH -t 20:00:00

module load launcher 
module load gsl

ulmit -S -m $((2*1024*1024))

export LD_LIBRARY_PATH=~/lib
export LAUNCHER_JOB_FILE=jobfile1
 
$LAUNCHER_DIR/paramrun


