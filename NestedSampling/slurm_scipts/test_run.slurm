

#SBATCH -J TestRun
#SBATCH -n 4
#SBATCH -N 1
#SBATCH -t 00:30:00

#create new ipython profile
profile=job_${SLURM_JOB_ID}
echo "Creating profile_${profile}"
ipython profile create ${profile}

#start controller
echo "Starting ipcontroller"
ipcontroller --ip="*" --profile=${profile} &
sleep 10

#run ipengines on each core
echo "Starting ipengines"
srun ipengine --profile=${profile} --location=$(hostname) &
sleep 25


