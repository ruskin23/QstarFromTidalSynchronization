#!/bin/bash

#SBATCH -J 6464285_emcee
#SBATCH -n 16 #Total number of tasks to run at one time
#SBATCH -N 1 #Number of nodes on which to distribute the tasks evenly
#SBATCH -o /home/rxp163130/scratch/version2/sampling_output/6464285.out
#SBATCH --mail-user=rxp163130@utdallas.edu
#SBATCH --mail-type=all
#SBATCH -t 4-00:00:00

export LD_LIBRARY_PATH=$LD_LIBRARAY_PATH:~/boost_lib:~/gsl_lib

export PYTHONPATH=/home/rxp163130/poet/PythonPackage:/home/rxp163130/poet/scripts

for d in /home/rxp163130/poet/build/libs/*/shared/release; do
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(readlink -e $d)
    export LIBRARY_PATH=$LIBRARY_PATH:$(readlink -e $d)
done


SOURCE_DIR=/home/rxp163130/QstarFromTidalSynchronization/MCMC/version2_emcee
OUTPUT_ROOT=/home/rxp163130/scratch/version2
OUTPUT_DIR=${OUTPUT_ROOT}/sampling_output


python3 \
    ${SOURCE_DIR}/mcmc_sampler.py \
    --system 6464285\
    --std-out-err-fname ${OUTPUT_DIR}'/%(system)s_%(now)s_%(pid)d.outerr'\
    --logging-fname ${OUTPUT_DIR}'/%(system)s_%(now)s_%(pid)d.log'\
    >\
    ${OUTPUT_DIR}/sample_6464285.outerr\
    2>&1
