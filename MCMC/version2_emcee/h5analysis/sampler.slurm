#!/bin/bash

#SBATCH -J sampler_param
#SBATCH -n 128  #Total number of tasks to run at one time
#SBATCH -N 1 #Number of nodes on which to distribute the tasks evenly
#SBATCH -o /work/06850/rpatel23/ls6/scratch/sampling_output/sampler/sampler.out
#SBATCH --mail-user=rxp163130@utdallas.edu
#SBATCH --mail-type=all
#SBATCH -t 00:30:00
#SBATCH -p development

module load intel/19.1.1
module load boost
module load gsl

export PYTHONPATH=/home1/06850/rpatel23/poet/PythonPackage:/home1/06850/rpatel23/poet/scripts

export PYTHONUNBUFFERED=1

for d in /home1/06850/rpatel23/poet/build/libs/*/shared/release; do
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(readlink -e $d)
    export LIBRARY_PATH=$LIBRARY_PATH:$(readlink -e $d)
done

SOURCE_DIR=/home1/06850/rpatel23/QstarFromTidalSynchronization/MCMC/version2_emcee/h5analysis
OUTPUT_ROOT=/work/06850/rpatel23/ls6/scratch
OUTPUT_DIR=${OUTPUT_ROOT}/sampling_output/sampler

python3 ${SOURCE_DIR}/samples_analysis.py --nprocs 128 --nsamples 2500 --npriors 1000 --pickle --sampler >${OUTPUT_DIR}/sampler.outerr 2>&1
