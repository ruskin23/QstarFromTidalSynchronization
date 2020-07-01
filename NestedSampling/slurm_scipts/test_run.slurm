#!/bin/bash

#SBATCH -J TestRun
#SBATCH -p debug
#SBATCH -n 4
#SBATCH -N 1
#SBATCH -t 00:30:00
#SBATCH -o /home/rxp163130/QstarFromTidalSynchronization/NestedSampling/slurm_scipts/output.out


python3 /work/06850/rpatel23/stampede2/QstarFromTidalSynchronization/NestedSampling/slurm_scipts/pathos_test_main.py -l 111
