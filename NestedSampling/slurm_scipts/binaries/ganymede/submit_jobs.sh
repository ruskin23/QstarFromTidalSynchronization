for S in 39,54,76,80,81,92,126
do
    sbatch  --job-name=nlive.$S --output=/home/rxp163130/scratch/Nested_Sampling/slurm_output.nlive_init.$S.ounlive_init.slurm --export=ALL,S=$S nlive_init.slurm
done