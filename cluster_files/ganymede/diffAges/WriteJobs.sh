#!/bin/bash

for  S in 2 5 14 18 23 27 33 35 37 49 59 62 63 65 68 69 77 78 82 99 100 111 117 119 125 128 129 132 133 134 138 142
do
	for  P in 10 20 30 40
		do
			echo "python3 /home/rxp163130/QstarFromTidalSynchronization/binary_star_evolution/analyze_spin_v_logQ/general_spin_v_logQ/diffages_main.py $S -p $P >/home/rxp163130/QstarFromTidalSynchronization/binary_star_evolution/analyze_spin_v_logQ/general_spin_v_logQ/output/ganymede/PercentileAges/"'$LAUNCHER_TSK_ID.out' 
	done
done >jobfile
