#!/bin/bash

#SBATCH -p normal
#SBATCH -J MCMCTest
#SBATCH -o /work/06850/rpatel23/stampede2/QstarFromTidalSynchronization/MCMC/combined/output/stampede2/test.out
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 00:15:00

#SBATCH --mail-user=rxp153130@utdallas.edu
#SBATCH --mail-type=all    

export PYTHONPATH=/work/06850/rpatel23/stampede2/poet/PythonPackage:/work/06850/rpatel23/stampede2/poet/scripts
export PYTHONUNBUFFERED=1

module load gsl
module load gcc/7.1.0
module load python3
module load boost/1.64

for d in /work/06850/rpatel23/stampede2/poet/build/libs/*/shared/release; do 
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(readlink -e $d)
    export LIBRARY_PATH=$LIBRARY_PATH:$(readlink -e $d)
done

python3 /work/06850/rpatel23/stampede2/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 1 -l 1 >/work/06850/rpatel23/stampede2/QstarFromTidalSynchronization/MCMC/combined/output/stampede2/test_1_1.txt 2>&1  
  
