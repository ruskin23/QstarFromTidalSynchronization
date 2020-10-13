#!/bin/bash



for I in {1..15}
do
    echo "./main.py -c -i $I -l 76 -m uncorrelated >/home/rxp163130/scratch/test_solver/76.$I.out"
done >jobfile
