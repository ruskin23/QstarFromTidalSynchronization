#!/bin/bash



for I in {11..20}
do
    echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -c -i $I -l 76 -m uncorrelated >/home/rxp163130/scratch/test_solver/76.$I.out"
done >jobfile
