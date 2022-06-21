#!/bin/bash

#SBATCH -J 4902030_emcee
#SBATCH -n 16 #Total number of tasks to run at one time
#SBATCH -N 1 #Number of nodes on which to distribute the tasks evenly
#SBATCH -o /home/rxp163130/scratch/version2/sampling_output/4902030.out
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
    --system 4902030\
    --std-out-err-fname ${OUTPUT_DIR}'/%(system)s_%(now)s_%(pid)d.outerr'\
    --logging-fname ${OUTPUT_DIR}'/%(system)s_%(now)s_%(pid)d.log'\
    >\
    ${OUTPUT_DIR}/sample_4902030.outerr\
    2>&1\
    &


# for SYSID in 10031409   10215422    10264202    10330495    10385682    10454401    10935310    10965961    1147276 11228612    11232745    11233911    11252617    11403216    11616200    11704044    12470530    2305543 3241344 3348093 3973504 4352168 4678171 4815612 4839180 5091614 5730394 7257373 7732791 7838639 7846730 8229048 8356054 8543278 8580438 8618226 8957954 10091257    10257903    10711551    10711913    1087921 10936427    10960995    10992733    11200773    11234677    11391181    12004679    12109845; do 
#     python3 \
#         ${SOURCE_DIR}/mcmc_sampler.py \
#         --system $SYSID\
#         --std-out-err-fname ${OUTPUT_DIR}'/%(system)s_%(now)s_%(pid)d.outerr'\
#         --logging-fname ${OUTPUT_DIR}'/%(system)s_%(now)s_%(pid)d.log'\
#         >\
#         ${OUTPUT_DIR}/sample_${SYSID}.outerr\
#         2>&1\
#         &
# done
