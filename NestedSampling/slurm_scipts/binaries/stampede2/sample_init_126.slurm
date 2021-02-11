#!/bin/bash

#SBATCH -p skx-normal
#SBATCH -J INIT126
#SBATCH -o /scratch/06850/rpatel23/Nested_Sampling/system_126/ slurm_output_11.out
#SBATCH -N 1
#SBATCH -n 44
#SBATCH -t 2-00:00:00

#SBATCH --mail-user=rxp153130@utdallas.edu
#SBATCH --mail-type=all

 

module load gsl
module load gcc
module load python3
module load boost

export PYTHONPATH=/work/06850/rpatel23/stampede2/poet/PythonPackage:/work/06850/rpatel23/stampede2/poet/scripts

for d in /work/06850/rpatel23/stampede2/poet/build/libs/*/shared/release; do
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(readlink -e $d)
    export LIBRARY_PATH=$LIBRARY_PATH:$(readlink -e $d)
done


python3 /work/06850/rpatel23/stampede2/QstarFromTidalSynchronization/NestedSampling/main/main.py -s continue -l 126 -n 45
