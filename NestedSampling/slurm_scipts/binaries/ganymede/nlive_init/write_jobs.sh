#!/bin/bash


MAIN_DIR="/home/rxp163130/QstarFromTidalSynchronization/NestedSampling/main"
OUT_DIR="/home/rxp163130/scratch/Nested_Sampling"

for S in 39 54 76 80 81 92 126
do
    for I in {1..32}
    do
        echo "python3 $MAIN_DIR/main.py -s start -l $S -i $I >$OUT_DIR/nlive_init.$S.$I.out"
    done >jobfile_$S
done 

