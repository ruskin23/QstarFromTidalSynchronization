#!/bin/bash


for S in 13 20 26 28 29 31 33 36 41 47 56 67 70 76 79 80 82 83 88 90 91 99 109 115 117 120 121 126 129 131 139 142
do
    echo "python3 /home/rxp163130/QstarFromTidalSynchronization/binary_star_evolution/analyze_spin_v_logQ/WindCatalog/main.py $S -b 0.0 -n >/home/rxp163130/scratch/bianry_evolution_analysis/$S.out" 
done >jobfile
