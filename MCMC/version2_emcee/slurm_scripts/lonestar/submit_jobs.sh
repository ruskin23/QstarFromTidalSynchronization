
for filename in *.slurm; do
sbatch $filename
sleep 5
done
