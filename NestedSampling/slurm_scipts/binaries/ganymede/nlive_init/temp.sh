#!/bin/bash

for S in 39 54 76 80 81 92 126
do
sbatch nlive_$S.slurm
done 