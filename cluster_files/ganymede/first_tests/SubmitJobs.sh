#!/bin/bash

for S in 1 8 12 13 25 32 36 39 
do
	for I in 1 2 3 4 5
	do
		sbatch --job-name=mcmc$S --output=output/output.$S.$I.out --export=S=$S,I=$I MCMC.sbatch 
	done
done
