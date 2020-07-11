#!/bin/bash




#SBATCH -J binaries 
#SBATCH -p debug
#SBATCH -n 1  #Total number of tasks to run at one time
#SBATCH -N 1  #Number of nodes on which to distribute the tasks evenly
#SBATCH -o /home/rxp163130/scratch/bianry_evolution_analysis/binaries.out
#SBATCH -t 00:30:00



module load gsl


ulimit -S -m $((2*1024*1024))

export LD_LIBRARY_PATH=~/lib



python3 /home/rxp163130/QstarFromTidalSynchronization/binary_star_evolution/analyze_spin_v_logQ/WindCatalog/main.py 13 -b 0.0 -n >/home/rxp163130/scratch/bianry_evolution_analysis/13.out

