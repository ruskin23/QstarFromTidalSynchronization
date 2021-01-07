#!/bin/bash

#SBATCH -J INIT80
#SBATCH -o /home/rxp163130/scratch/Nested_Sampling/system_80/slurm_output_1.out
#SBATCH -N 3
#SBATCH -n 45
#SBATCH -t 4-00:00:00
 
ulimit -S -m $((2*1024*1024))

module load launcher

export LD_LIBRARY_PATH=$LD_LIBRARAY_PATH:~/boost_lib:~/gsl_lib

export PYTHONPATH=/home/rxp163130/poet/PythonPackage:/home/rxp163130/poet/scripts

for d in /home/rxp163130/poet/build/libs/*/shared/release; do
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(readlink -e $d)
    export LIBRARY_PATH=$LIBRARY_PATH:$(readlink -e $d)
done

python3 /home/rxp163130/QstarFromTidalSynchronization/NestedSampling/main/main.py -s continue -l 80 -n 45
