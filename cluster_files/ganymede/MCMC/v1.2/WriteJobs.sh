#!/bin/bash
#low acceptance ratio:
#39 137 54 80 126 76 1 32 67 81 95 96 50 85 56 73 86 92 20
#remaining:
#17 25 109 120 83 106 57 36 88 44 123 28 47 79 94 13 48 43 93 70 8 12 31
#diff:
#111 129 132 18 27 33 49 5 62 68 82 99 

rm uncorrelated

for B in 1 2 3 4 5
do
    echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -c -i $B -l 1 -m adaptive >/home/rxp163130/scratch/tenth_run/1.$B.out"
	for A in 39 137 54 80 126 76 1 32 67 81 95 96 50 85 56 73 86 92 20
	do
		echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -c -i $B -l $A -m adaptive >/home/rxp163130/scratch/tenth_run/$A.$B.out"
	done
done >adaptive



for B in 1 2 3 4 5
do
    for A in 17 25 120 83 57 36 88 44 47 79 94 13 48 43 93 70 8 31
    do
        echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -c -i $B -l $A -m uncorrelated >/home/rxp163130/scratch/tenth_run/$A.$B.out"
    done
done >uncorrelated

#for B in 1 2 3 4 5
#do
#    for A in 111 128 132 18 27 33 49 5 62 68 99 
#    do
#        echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -c -i $B -l $A -p $B -m uncorrelated >/home/rxp163130/scratch/tenth_run/$A.$B.out" >> uncorrelated
#    done
#done 



