#!/bin/bash


for B in 1.0 2.0 3.0 4.0
do
	for S in 1 8 12 13
	do
		echo "python3 /home/rxp163130/QstarFromTidalSynchronization/binary_star_evolution/updated_evolution_code/ganymede_main.py -s $S -b $B >/home/rxp163130/scratch/bianry_evolution_analysis/$S.$B.out"
	done
done >jobfile1

for B in 1.0 2.0 3.0 4.0
do
	for S in 17 20 25 28 
	do
		echo "python3 /home/rxp163130/QstarFromTidalSynchronization/binary_star_evolution/updated_evolution_code/ganymede_main.py -s $S -b $B >/home/rxp163130/scratch/bianry_evolution_analysis/$S.$B.out"
	done
done >jobfile2

for B in 1.0 2.0 3.0 4.0
do
	for S in 32 36 39 43
	do
		echo "python3 /home/rxp163130/QstarFromTidalSynchronization/binary_star_evolution/updated_evolution_code/ganymede_main.py -s $S -b $B >/home/rxp163130/scratch/bianry_evolution_analysis/$S.$B.out"
	done
done >jobfile3

for B in 1.0 2.0 3.0 4.0
do
	for S in 44 47 48 50
	do
		echo "python3 /home/rxp163130/QstarFromTidalSynchronization/binary_star_evolution/updated_evolution_code/ganymede_main.py -s $S -b $B >/home/rxp163130/scratch/bianry_evolution_analysis/$S.$B.out"
	done
done >jobfile4

for B in 1.0 2.0 3.0 4.0
do
	for S in 54 56 67 70
	do
		echo "python3 /home/rxp163130/QstarFromTidalSynchronization/binary_star_evolution/updated_evolution_code/ganymede_main.py -s $S -b $B >/home/rxp163130/scratch/bianry_evolution_analysis/$S.$B.out"
	done
done >jobfile5

for B in 1.0 2.0 3.0 4.0
do
	for S in 73 76 79 80
	do
		echo "python3 /home/rxp163130/QstarFromTidalSynchronization/binary_star_evolution/updated_evolution_code/ganymede_main.py -s $S -b $B >/home/rxp163130/scratch/bianry_evolution_analysis/$S.$B.out"
	done
done >jobfile6

for B in 1.0 2.0 3.0 4.0
do
	for S in 81 83 84 85
	do
		echo "python3 /home/rxp163130/QstarFromTidalSynchronization/binary_star_evolution/updated_evolution_code/ganymede_main.py -s $S -b $B >/home/rxp163130/scratch/bianry_evolution_analysis/$S.$B.out"
	done
done >jobfile7

for B in 1.0 2.0 3.0 4.0
do
	for S in 86 88 92 93
	do
		echo "python3 /home/rxp163130/QstarFromTidalSynchronization/binary_star_evolution/updated_evolution_code/ganymede_main.py -s $S -b $B >/home/rxp163130/scratch/bianry_evolution_analysis/$S.$B.out"
	done
done >jobfile8

for B in 1.0 2.0 3.0 4.0
do
	for S in 94 95 96 106
	do
		echo "python3 /home/rxp163130/QstarFromTidalSynchronization/binary_star_evolution/updated_evolution_code/ganymede_main.py -s $S -b $B >/home/rxp163130/scratch/bianry_evolution_analysis/$S.$B.out"
	done
done >jobfile9

for B in 1.0 2.0 3.0 4.0
do
	for S in 109 120 123 126
	do
		echo "python3 /home/rxp163130/QstarFromTidalSynchronization/binary_star_evolution/updated_evolution_code/ganymede_main.py -s $S -b $B >/home/rxp163130/scratch/bianry_evolution_analysis/$S.$B.out"
	done
done >jobfile10

for B in 1.0 2.0 3.0 4.0
do
	for S in 137
	do
		echo "python3 /home/rxp163130/QstarFromTidalSynchronization/binary_star_evolution/updated_evolution_code/ganymede_main.py -s $S -b $B >/home/rxp163130/scratch/bianry_evolution_analysis/$S.$B.out"
	done
done >jobfile11
