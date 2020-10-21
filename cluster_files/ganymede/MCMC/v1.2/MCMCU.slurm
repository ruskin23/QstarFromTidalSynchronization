#!/bin/bash

#SBATCH -J MCMCU 
#SBATCH -n 65  #Total number of tasks to run at one time
#SBATCH -N 5  #Number of nodes on which to distribute the tasks evenly
#SBATCH -o /home/rxp163130/scratch/Run28/MCMCU.out
#SBATCH -t 4-00:00:00


ulimit -S -m $((2*1024*1024))

module load launcher

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/boost_lib:~/gsl_lib
export PYTHONPATH=/home/rxp163130/poet/PythonPackage:/home/rxp163130/poet/scripts

for d in /home/rxp163130/poet/build/libs/*/shared/release; do
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(readlink -e $d)
    export LIBRARY_PATH=$LIBRARY_PATH:$(readlink -e $d)
done

export LAUNCHER_JOB_FILE=uncorrelated 
$LAUNCHER_DIR/paramrun


