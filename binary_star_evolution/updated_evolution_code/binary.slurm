#!/bin/bash

#SBATCH -J binary_evolution
#SBATCH -n 27  #Total number of tasks to run at one time
#SBATCH -N 1 #Number of nodes on which to distribute the tasks evenly
#SBATCH -o /work/06850/rpatel23/ls6/scratch/binary.out
#SBATCH --mail-user=rxp163130@utdallas.edu
#SBATCH --mail-type=all
#SBATCH -t 02:00:00
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

python3 Calculate_Final_Spin.py >/work/06850/rpatel23/ls6/scratch/binary_root_hybr.outerr -f root -m hybr &
python3 Calculate_Final_Spin.py >/work/06850/rpatel23/ls6/scratch/binary_root_lm.outerr -f root -m lm &
python3 Calculate_Final_Spin.py >/work/06850/rpatel23/ls6/scratch/binary_root_broyden1.outerr -f root -m broyden1 &
python3 Calculate_Final_Spin.py >/work/06850/rpatel23/ls6/scratch/binary_root_broyden2.outerr -f root -m broyden2 &
python3 Calculate_Final_Spin.py >/work/06850/rpatel23/ls6/scratch/binary_root_anderson.outerr -f root -m anderson &
python3 Calculate_Final_Spin.py >/work/06850/rpatel23/ls6/scratch/binary_root_linearmixing.outerr -f root -m linearmixing &
python3 Calculate_Final_Spin.py >/work/06850/rpatel23/ls6/scratch/binary_root_diagbroyden.outerr -f root -m diagbroyden &
python3 Calculate_Final_Spin.py >/work/06850/rpatel23/ls6/scratch/binary_root_excitingmixing.outerr -f root -m excitingmixing &
python3 Calculate_Final_Spin.py >/work/06850/rpatel23/ls6/scratch/binary_root_krylov.outerr -f root -m krylov &
python3 Calculate_Final_Spin.py >/work/06850/rpatel23/ls6/scratch/binary_root_df.outerr -f root -m df-sane &
python3 Calculate_Final_Spin.py >/work/06850/rpatel23/ls6/scratch/binary_minimize_Nelder-Mead.outerr -f minimize -m Nelder-Mead &
python3 Calculate_Final_Spin.py >/work/06850/rpatel23/ls6/scratch/binary_minimize_Powell.outerr -f minimize -m Powell &
python3 Calculate_Final_Spin.py >/work/06850/rpatel23/ls6/scratch/binary_minimize_CG.outerr -f minimize -m CG &
python3 Calculate_Final_Spin.py >/work/06850/rpatel23/ls6/scratch/binary_minimize_BFGS.outerr -f minimize -m BFGS &
python3 Calculate_Final_Spin.py >/work/06850/rpatel23/ls6/scratch/binary_minimize_Newton-CG.outerr -f minimize -m Newton-CG &
python3 Calculate_Final_Spin.py >/work/06850/rpatel23/ls6/scratch/binary_minimize_L-BFGS-B.outerr -f minimize -m L-BFGS-B &
python3 Calculate_Final_Spin.py >/work/06850/rpatel23/ls6/scratch/binary_minimize__minimize_TNC.outerr -f minimize -m TNC &
python3 Calculate_Final_Spin.py >/work/06850/rpatel23/ls6/scratch/binary_minimize_COBYLA.outerr -f minimize -m COBYLA &
python3 Calculate_Final_Spin.py >/work/06850/rpatel23/ls6/scratch/binary_minimize_SLSQP.outerr -f minimize -m SLSQP &
python3 Calculate_Final_Spin.py >/work/06850/rpatel23/ls6/scratch/binary_minimize_trust-constr.outerr -f minimize -m trust-constr &
python3 Calculate_Final_Spin.py >/work/06850/rpatel23/ls6/scratch/binary_minimize_dogleg.outerr -f minimize -m dogleg &
python3 Calculate_Final_Spin.py >/work/06850/rpatel23/ls6/scratch/binary_minimize_trust-ncg.outerr -f minimize -m trust-ncg &
python3 Calculate_Final_Spin.py >/work/06850/rpatel23/ls6/scratch/binary_minimize_trust-exact.outerr -f minimize -m trust-exact &
python3 Calculate_Final_Spin.py >/work/06850/rpatel23/ls6/scratch/binary_minimize_trust-krylov.outerr -f minimize -m trust-krylov &
python3 Calculate_Final_Spin.py >/work/06850/rpatel23/ls6/scratch/binary_least_squares_trf.outerr -f least_squares -m trf &
python3 Calculate_Final_Spin.py >/work/06850/rpatel23/ls6/scratch/binary_least_squares_dogbox.outerr -f least_squares -m dogbox &
python3 Calculate_Final_Spin.py >/work/06850/rpatel23/ls6/scratch/binary_least_squares_lm.outerr -f least_squares -m lm &


wait