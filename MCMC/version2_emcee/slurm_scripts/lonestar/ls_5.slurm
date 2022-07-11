#!/bin/bash

#SBATCH -J emcee_5
#SBATCH -n 128  #Total number of tasks to run at one time
#SBATCH -N 1 #Number of nodes on which to distribute the tasks evenly
#SBATCH -o /work/06850/rpatel23/ls6/scratch/sampling_output/emcee_5.out
#SBATCH --mail-user=rxp163130@utdallas.edu
#SBATCH --mail-type=all
#SBATCH -t 2-00:00:00
#SBATCH -p normal

module load intel/19.1.1
module load boost
module load gsl

export PYTHONPATH=/home1/06850/rpatel23/poet/PythonPackage:/home1/06850/rpatel23/poet/scripts

export PYTHONUNBUFFERED=1

run_system() {
    SOURCE_DIR=/home1/06850/rpatel23/QstarFromTidalSynchronization/MCMC/version2_emcee
    OUTPUT_ROOT=/work/06850/rpatel23/ls6/scratch
    OUTPUT_DIR=${OUTPUT_ROOT}/sampling_output
    sysid=$1

    while true; do
        python3 \
            ${SOURCE_DIR}/mcmc_sampler.py \
            --system ${sysid}\
            --std_out_err_fname ${OUTPUT_DIR}'/%(system)s_%(now)s_%(pid)d.outerr'\
            --logging_fname ${OUTPUT_DIR}'/%(system)s_%(now)s_%(pid)d.log'\
            >\
            ${OUTPUT_DIR}/sample_${sysid}.outerr\
            2>&1
        done
}

for SYSID in 10960995   10992733    11200773    11391181   12004679 12351927    2447893 2449090; do 
    run_system $SYSID & 
done

wait
