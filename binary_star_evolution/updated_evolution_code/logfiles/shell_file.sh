#!/bin/bash

if [ -e parsed.csv ]
then
    rm parsed.csv
fi

for I in 29;do
    echo $I
    tar -xf newRun$I.tar.gz
    for S in system_*; do
        python3 data_processing.py $S    
    done
    rm -rf system_*
done

