
for filename in *.slurm; do
sbatch $filename
wait 1
done
