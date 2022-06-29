#!/bin/bash

#SBATCH -J 5288543_emcee
#SBATCH -n 16 #Total number of tasks to run at one time
#SBATCH -N 1 #Number of nodes on which to distribute the tasks evenly
#SBATCH -o /home/rxp163130/scratch/version2/sampling_output/5288543.out
#SBATCH --mail-user=rxp163130@utdallas.edu
#SBATCH --mail-type=all
#SBATCH -t 4-00:00:00

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/rxp163130/lib
export LIBRARY_PATH=$LIBRARY_PATH:/home/rxp163130/lib


export PYTHONPATH=/home/rxp163130/poet/PythonPackage:/home/rxp163130/poet/scripts

export PYTHONUNBUFFERED=1

SOURCE_DIR=/home/rxp163130/QstarFromTidalSynchronization/MCMC/version2_emcee
OUTPUT_ROOT=/home/rxp163130/scratch/version2
OUTPUT_DIR=${OUTPUT_ROOT}/sampling_output


python3 \
    ${SOURCE_DIR}/mcmc_sampler.py \
    --system 5288543\
    --std_out_err_fname ${OUTPUT_DIR}'/%(system)s_%(now)s_%(pid)d.outerr'\
    --logging_fname ${OUTPUT_DIR}'/%(system)s_%(now)s_%(pid)d.log'\
    >\
    ${OUTPUT_DIR}/sample_5288543.outerr\
    2>&1
