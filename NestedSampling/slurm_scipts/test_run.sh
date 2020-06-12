#!/bin/bash



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

