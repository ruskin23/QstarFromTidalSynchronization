#!/bin/bash


for B in 1 2 3 4 5
do
    echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -c -i $B -l 1 -m adaptive >/home/rxp163130/scratch/Run20/1.$B.out"
	for A in 39 137 54 80 126 76 1 32 67 81 95 96 85 56 73 86 92 20
	do
		echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -c -i $B -l $A -m adaptive >/home/rxp163130/scratch/Run20/$A.$B.out"
	done
done >adaptive



for B in 1 2 3 4 5
do
    for A in 17 83 57 36 88 44 94 13 48 43 93 70 31
    do
        echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -c -i $B -l $A -m uncorrelated >/home/rxp163130/scratch/Run20/$A.$B.out"
    done
done >uncorrelated

