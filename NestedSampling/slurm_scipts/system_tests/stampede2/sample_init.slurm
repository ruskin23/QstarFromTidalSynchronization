#!/bin/bash

#SBATCH -p skx-normal
#SBATCH -J Nested84
#SBATCH -n 90
#SBATCH -N 2
#SBATCH -t 00:30:00
#SBATCH -o /scratch/06850/rpatel23/output/NestedSampling/test_sampleinit.out

module load gsl
module load python3

ulimit -S -m $((2*1024*1024))

export LD_LIBRARY_PATH=~/lib

python3 /work/06850/rpatel23/stampede2/QstarFromTidalSynchronization/NestedSampling/slurm_scipts/system_tests/main.py -l 84 -t Internally >/scratch/06850/rpatel23/output/NestedSampling/84_sampleinit_test.txt
