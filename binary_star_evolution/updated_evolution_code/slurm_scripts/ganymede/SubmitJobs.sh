#!/bin/bash


for J in 1 2 3 4 5 6 7 8 9 10 11
do
	sbatch binary_scipts_$J.slurm
	sleep 10
done
