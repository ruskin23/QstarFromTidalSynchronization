#!/bin/bash


for B in -3.0 -2.0 -1.0 1.0 2.0 3.0
do
	for S in 1 8 12 13 17 20 25 28 32 36 39 43 44 47 48 50 54 56 67 70 73 76 79 80 81 83 84 85 86 88 92 93 94 95 96 106 109 120 123 126 137
	do
		echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/version1/ganymede_main.py -s $S -b $B >/home/rxp163130/scratch/bianry_evolution_analysis/$S.$B.out"
	done
done >jobfile


