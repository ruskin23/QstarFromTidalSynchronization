#!/bin/bash

#SBATCH -J emcee_12
#SBATCH -n 128  #Total number of tasks to run at one time
#SBATCH -N 8 #Number of nodes on which to distribute the tasks evenly
#SBATCH -o /work/06850/rpatel23/ls6/scratch/sampling_output/emcee_12.out
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

for SYS in 6778050 6863840 6927629 6949550 6962018 7025540 7336754 7362852;
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
