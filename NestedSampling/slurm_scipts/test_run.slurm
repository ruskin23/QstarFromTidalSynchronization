#!/bin/bash

#SBATCH -J TestRun
#SBATCH -p debug
#SBATCH -n 4
#SBATCH -N 1
#SBATCH -t 00:30:00
#SBATCH -o /home/rxp163130/QstarFromTidalSynchronization/NestedSampling/slurm_scipts/output.out 

profile=job_${SLURM_JOB_ID}_$(hostname)

echo "Creating profile ${profile}"
ipython profile create ${profile}

echo "Launching controller"
ipcontroller --ip="*" --profile=${profile} --log-to-file &
sleep 5

echo "Launching engines"
srun ipengine --profile=${profile} --location=$(hostname) --log-to-file &
sleep 5


python3 test_main.py -l 1 -p ${profile}
