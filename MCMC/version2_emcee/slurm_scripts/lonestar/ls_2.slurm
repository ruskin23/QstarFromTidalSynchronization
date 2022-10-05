#!/bin/bash

#SBATCH -J emcee_2
#SBATCH -n 128  #Total number of tasks to run at one time
#SBATCH -N 1 #Number of nodes on which to distribute the tasks evenly
#SBATCH -o /work/06850/rpatel23/ls6/scratch/sampling_output/node_2/emcee_2.out
#SBATCH --mail-user=rxp163130@utdallas.edu
#SBATCH --mail-type=all
#SBATCH -t 2-00:00:00
#SBATCH -p normal

module load intel/19.1.1
module load boost
module load gsl

export PYTHONPATH=/home1/06850/rpatel23/poet/PythonPackage:/home1/06850/rpatel23/poet/scripts

export PYTHONUNBUFFERED=1

for d in /home1/06850/rpatel23/poet/build/libs/*/shared/release; do
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(readlink -e $d)
    export LIBRARY_PATH=$LIBRARY_PATH:$(readlink -e $d)
done

run_system() {
    sysid=$1
    SOURCE_DIR=/home1/06850/rpatel23/QstarFromTidalSynchronization/MCMC/version2_emcee
    OUTPUT_ROOT=/work/06850/rpatel23/ls6/scratch
    OUTPUT_DIR=${OUTPUT_ROOT}/sampling_output/node_2/system_${sysid}

    if [ ! -d ${OUTPUT_DIR} ] 
    then
        mkdir ${OUTPUT_DIR}
    fi

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



for SYSID in 8618226 8957954 10992733 4252226 4285087 4579321 10385682 4380283; do 
    run_system $SYSID & 
done

wait
