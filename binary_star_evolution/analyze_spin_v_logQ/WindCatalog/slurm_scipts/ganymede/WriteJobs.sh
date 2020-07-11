#!/bin/bash


for S in 1 2 5 8 10 12 13 14 17 20 25 26 27 28 29 31 32 33 36 37 41 43 44 47 48 49 52 54 55 56 57 59 61 62 63 67 69 70 73 76 77 78 79 80 81 82 83 84 88 90 91 93 95 99 100 101 103 105 106 107 109 110 112 113 114 115 117 120 121 124 125 126 128 129 131 133 134 137 138 139 142
do
    echo "python3 /home/rxp163130/QstarFromTidalSynchronization/binary_star_evolution/analyze_spin_v_logQ/WindCatalog/main.py $S -b 0.0 -n >/home/rxp163130/scratch/bianry_evolution_analysis/$S.out" 
done >jobfile
