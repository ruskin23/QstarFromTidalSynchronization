#!/bin/bash

#SBATCH -J emcee_2
#SBATCH -n 128  #Total number of tasks to run at one time
#SBATCH -N 8 #Number of nodes on which to distribute the tasks evenly
#SBATCH -o /work/06850/rpatel23/ls6/scratch/sampling_output/emcee_2.out
#SBATCH --mail-user=rxp163130@utdallas.edu
#SBATCH --mail-type=all
#SBATCH -t 2-00:00:00

module laod boost
module load gsl

export PYTHONPATH=/home1/06850/rpatel23/poet/PythonPackage:/home1/06850/rpatel23/poet/scripts

export PYTHONUNBUFFERED=1

SOURCE_DIR=/home1/06850/rpatel23/QstarFromTidalSynchronization/MCMC/version2_emcee
OUTPUT_ROOT=/work/06850/rpatel23/ls6/scratch
OUTPUT_DIR=${OUTPUT_ROOT}/sampling_output

for SYS in 12470530 2305543 3241344 3348093 3973504 4352168 4678171 4815612;
do 
python3 \
    ${SOURCE_DIR}/mcmc_sampler.py \
    --system ${SYS}\
    --std_out_err_fname ${OUTPUT_DIR}'/%(system)s_%(now)s_%(pid)d.outerr'\
    --logging_fname ${OUTPUT_DIR}'/%(system)s_%(now)s_%(pid)d.log'\
    >\
    ${OUTPUT_DIR}/sample_${SYS}.outerr\
    2>&1\
    &
done
