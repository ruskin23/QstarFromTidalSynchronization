#!/bin/bash

#SBATCH -p normal
#SBATCH -J Interpolator
#SBATCH -o interpolation_output.out
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 16:00:00

#SBATCH --mail-user=rxp153130@utdallas.edu
#SBATCH --mail-type=all    

export PYTHONPATH=/work/06850/rpatel23/stampede2/poet/PythonPackage:/work/06850/rpatel23/stampede2/poet/scripts

module load gsl
module load gcc
module load python3
module load boost

for d in /work/06850/rpatel23/stampede2/poet/build/libs/*/shared/release; do 
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(readlink -e $d)
    export LIBRARY_PATH=$LIBRARY_PATH:$(readlink -e $d)
done

python3 /work/06850/rpatel23/stampede2/poet/scripts/add_default_interpolator.py -r '/work/06850/rpatel23/stampede2/poet/stellar_interpolators' -t '/work/06850/rpatel23/stampede2/poet/MESA_tracks' --threads 1
  
