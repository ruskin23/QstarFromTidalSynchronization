#!/bin/bash

#SBATCH -p skx-normal
#SBATCH -J DNested84
#SBATCH -n 45
#SBATCH -N 1
#SBATCH -t 2-00:00:00
#SBATCH -o /scratch/06850/rpatel23/output/NestedSampling/84.dynamicCompleteSlurm.out

module load gsl/2.3
module load gcc
module load python3
module load boost


ulimit -S -m $((2*1024*1024))

export PYTHONPATH=/work/06850/rpatel23/stampede2/poet/PythonPackage:/work/06850/rpatel23/stampede2/poet/scripts

for d in /work/06850/rpatel23/stampede2/poet/build/libs/*/shared/release; do
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(readlink -e $d)
    export LIBRARY_PATH=$LIBRARY_PATH:$(readlink -e $d)
done

python3 /work/06850/rpatel23/stampede2/QstarFromTidalSynchronization/NestedSampling/slurm_scipts/system_tests/main.py -l 84 -t Internally >/scratch/06850/rpatel23/output/NestedSampling/84_dynamic_complete.txt
