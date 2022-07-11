#!/bin/bash

#SBATCH -J emcee_6
#SBATCH -n 128  #Total number of tasks to run at one time
#SBATCH -N 8 #Number of nodes on which to distribute the tasks evenly
#SBATCH -o /work/06850/rpatel23/ls6/scratch/sampling_output/emcee_6.out
#SBATCH --mail-user=rxp163130@utdallas.edu
#SBATCH --mail-type=all
#SBATCH -t 2-00:00:00

module load intel/19.1.1
module load boost
module load gsl

export PYTHONPATH=/home1/06850/rpatel23/poet/PythonPackage:/home1/06850/rpatel23/poet/scripts

export PYTHONUNBUFFERED=1


run_system() {
    SOURCE_DIR=/home1/06850/rpatel23/QstarFromTidalSynchronization/MCMC/version2_emcee
    OUTPUT_ROOT=/work/06850/rpatel23/ls6/scratch
    OUTPUT_DIR=${OUTPUT_ROOT}/sampling_output

    while true; do
        python3 \
            ${SOURCE_DIR}/mcmc_sampler.py \
            --system ${SYS}\
            --std_out_err_fname ${OUTPUT_DIR}'/%(system)s_%(now)s_%(pid)d.outerr'\
            --logging_fname ${OUTPUT_DIR}'/%(system)s_%(now)s_%(pid)d.log'\
            >\
            ${OUTPUT_DIR}/sample_${SYS}.outerr\
            2>&1
        done
}

for SYSID in 2852560 2860788 3098194 3248019 3323289 3344427 3834364 3838496; do 
    run_system $SYSID & 
done

wait