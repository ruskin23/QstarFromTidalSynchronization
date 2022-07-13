#!/bin/bash

#SBATCH -J binary_evolution
#SBATCH -n 1  #Total number of tasks to run at one time
#SBATCH -N 1 #Number of nodes on which to distribute the tasks evenly
#SBATCH -o /work/06850/rpatel23/ls6/scratch/binary.out
#SBATCH --mail-user=rxp163130@utdallas.edu
#SBATCH --mail-type=all
#SBATCH -t 02:00:00
#SBATCH -p debug

module load intel/19.1.1
module load boost
module load gsl

export PYTHONPATH=/home1/06850/rpatel23/poet/PythonPackage:/home1/06850/rpatel23/poet/scripts

export PYTHONUNBUFFERED=1

for d in /home1/06850/rpatel23/poet/build/libs/*/shared/release; do
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(readlink -e $d)
    export LIBRARY_PATH=$LIBRARY_PATH:$(readlink -e $d)
done

python3 Calculate_Final_Spin.py >/work/06850/rpatel23/ls6/scratch/binary.outerr
