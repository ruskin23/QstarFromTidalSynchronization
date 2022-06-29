#!/bin/bash

#SBATCH -J 7374746_emcee
#SBATCH -n 16 #Total number of tasks to run at one time
#SBATCH -N 1 #Number of nodes on which to distribute the tasks evenly
#SBATCH -o /work/06850/rpatel23/ls6/scratch/sampling_output/debug/7374746.out
#SBATCH --mail-user=rxp163130@utdallas.edu
#SBATCH --mail-type=all
#SBATCH -t 2:00:00
#SBATCH -p debug

module laod boost
module load gsl

export PYTHONPATH=/home1/06850/rpatel23/poet/PythonPackage:/home1/06850/rpatel23/poet/scripts

export PYTHONUNBUFFERED=1

SOURCE_DIR=/home1/06850/rpatel23/QstarFromTidalSynchronization/MCMC/version2_emcee
OUTPUT_ROOT=/work/06850/rpatel23/ls6/scratch
OUTPUT_DIR=${OUTPUT_ROOT}/sampling_output/debug


python3 \
    ${SOURCE_DIR}/mcmc_sampler.py \
    --system 7374746\
    --std_out_err_fname ${OUTPUT_DIR}'/%(system)s_%(now)s_%(pid)d.outerr'\
    --logging_fname ${OUTPUT_DIR}'/%(system)s_%(now)s_%(pid)d.log'\
    >\
    ${OUTPUT_DIR}/sample_7374746.outerr\
    2>&1
