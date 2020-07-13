#!/bin/bash

for B in 1 2 3 4 5
do
	for A in 8 12 13 20 26 28 29 31 33 36 41 43 47 56 63 67 70 76 79 80 82 83 88 90 91 99 109 115 117 120 121 126 128 129 131 139 142
	do
		echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/version2/main.py -s -i $B -l $A >/home/rxp163130/scratch/version2/first_run/$A.$B.out"
	done
done >jobfile


