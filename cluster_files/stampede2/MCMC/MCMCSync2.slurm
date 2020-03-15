#!/bin/bash

#SBATCH -p skx-normal
#SBATCH -J MCMCNonSync2
#SBATCH -o /scratch/06850/rpatel23/output/MCMCS2.o%j
#SBATCH -N 1
#SBATCH -n 45
#SBATCH -t 2-00:00:00

#SBATCH --mail-user=rxp153130@utdallas.edu
#SBATCH --mail-type=all    

ulimit -S -m $((2*1024*1024))

module load gsl
module load gcc
module load python3
module load boost
module load launcher

export PYTHONPATH=/work/06850/rpatel23/stampede2/poet/PythonPackage:/work/06850/rpatel23/stampede2/poet/scripts

for d in /work/06850/rpatel23/stampede2/poet/build/libs/*/shared/release; do 
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(readlink -e $d)
    export LIBRARY_PATH=$LIBRARY_PATH:$(readlink -e $d)
done

export LAUNCHER_JOB_FILE=SyncJob2

$LAUNCHER_DIR/paramrun
