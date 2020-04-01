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
    echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -c -i $B -l 1 -m adaptive >/home/rxp163130/scratch/fifth_run/1.$B.out"
	for A in 39 137 54 80 126 76 1 32 67 81 95 96 50 85 56 73 86 92 20
	do
		echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -c -i $B -l $A -m adaptive >/home/rxp163130/scratch/fifth_run/$A.$B.out"
	done
done >adaptive



for B in 1 2 3 4 5
do
    for A in 17 25 109 120 83 106 57 36 88 44 123 28 47 79 94 13 48 43 93 70 8 12 31
    do
        echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -c -i $B -l $A -m uncorrelated >/home/rxp163130/scratch/fifth_run/$A.$B.out"
    done
done >uncorrelated

echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 1 -l 111 -p 2 -m uncorrelated >/home/rxp163130/scratch/fifth_run/111.1.out" >> uncorrelated
echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 2 -l 111 -p 2 -m uncorrelated >/home/rxp163130/scratch/fifth_run/111.2.out" >> uncorrelated
echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 3 -l 111 -p 2 -m uncorrelated >/home/rxp163130/scratch/fifth_run/111.3.out" >> uncorrelated
echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 4 -l 111 -p 3 -m uncorrelated >/home/rxp163130/scratch/fifth_run/111.4.out" >> uncorrelated
echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 5 -l 111 -p 3 -m uncorrelated >/home/rxp163130/scratch/fifth_run/111.5.out" >> uncorrelated

echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 1 -l 129 -p 1 -m uncorrelated >/home/rxp163130/scratch/fifth_run/129.1.out" >> uncorrelated
echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 2 -l 129 -p 2 -m uncorrelated >/home/rxp163130/scratch/fifth_run/129.2.out" >> uncorrelated
echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 3 -l 129 -p 2 -m uncorrelated >/home/rxp163130/scratch/fifth_run/129.3.out" >> uncorrelated
echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 4 -l 129 -p 3 -m uncorrelated >/home/rxp163130/scratch/fifth_run/129.4.out" >> uncorrelated
echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 5 -l 129 -p 3 -m uncorrelated >/home/rxp163130/scratch/fifth_run/129.5.out" >> uncorrelated

echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 1 -l 132 -p 3 -m uncorrelated >/home/rxp163130/scratch/fifth_run/132.1.out" >> uncorrelated
echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 2 -l 132 -p 3 -m uncorrelated >/home/rxp163130/scratch/fifth_run/132.2.out" >> uncorrelated
echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 3 -l 132 -p 4 -m uncorrelated >/home/rxp163130/scratch/fifth_run/132.3.out" >> uncorrelated
echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 4 -l 132 -p 4 -m uncorrelated >/home/rxp163130/scratch/fifth_run/132.4.out" >> uncorrelated
echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 5 -l 132 -p 4 -m uncorrelated >/home/rxp163130/scratch/fifth_run/132.5.out" >> uncorrelated


echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 1 -l 18 -p 2 -m uncorrelated >/home/rxp163130/scratch/fifth_run/18.1.out" >> uncorrelated
echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 2 -l 18 -p 3 -m uncorrelated >/home/rxp163130/scratch/fifth_run/18.2.out" >> uncorrelated
echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 3 -l 18 -p 3 -m uncorrelated >/home/rxp163130/scratch/fifth_run/18.3.out" >> uncorrelated
echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 4 -l 18 -p 4 -m uncorrelated >/home/rxp163130/scratch/fifth_run/18.4.out" >> uncorrelated
echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 5 -l 18 -p 4 -m uncorrelated >/home/rxp163130/scratch/fifth_run/18.5.out" >> uncorrelated

echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 1 -l 27 -p 1 -m uncorrelated >/home/rxp163130/scratch/fifth_run/27.1.out" >> uncorrelated
echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 2 -l 27 -p 1 -m uncorrelated >/home/rxp163130/scratch/fifth_run/27.2.out" >> uncorrelated
echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 3 -l 27 -p 1 -m uncorrelated >/home/rxp163130/scratch/fifth_run/27.3.out" >> uncorrelated
echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 4 -l 27 -p 1 -m uncorrelated >/home/rxp163130/scratch/fifth_run/27.4.out" >> uncorrelated
echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 5 -l 27 -p 1 -m uncorrelated >/home/rxp163130/scratch/fifth_run/27.5.out" >> uncorrelated

echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 1 -l 33 -p 2 -m uncorrelated >/home/rxp163130/scratch/fifth_run/33.1.out" >> uncorrelated
echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 2 -l 33 -p 3 -m uncorrelated >/home/rxp163130/scratch/fifth_run/33.2.out" >> uncorrelated
echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 3 -l 33 -p 3 -m uncorrelated >/home/rxp163130/scratch/fifth_run/33.3.out" >> uncorrelated
echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 4 -l 33 -p 4 -m uncorrelated >/home/rxp163130/scratch/fifth_run/33.4.out" >> uncorrelated
echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 5 -l 33 -p 4 -m uncorrelated >/home/rxp163130/scratch/fifth_run/33.5.out" >> uncorrelated

echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 1 -l 49 -p 2 -m uncorrelated >/home/rxp163130/scratch/fifth_run/49.1.out" >> uncorrelated
echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 2 -l 49 -p 3 -m uncorrelated >/home/rxp163130/scratch/fifth_run/49.2.out" >> uncorrelated
echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 3 -l 49 -p 3 -m uncorrelated >/home/rxp163130/scratch/fifth_run/49.3.out" >> uncorrelated
echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 4 -l 49 -p 4 -m uncorrelated >/home/rxp163130/scratch/fifth_run/49.4.out" >> uncorrelated
echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 5 -l 49 -p 4 -m uncorrelated >/home/rxp163130/scratch/fifth_run/49.5.out" >> uncorrelated

echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 1 -l 5 -p 3 -m uncorrelated >/home/rxp163130/scratch/fifth_run/5.1.out" >> uncorrelated
echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 2 -l 5 -p 4 -m uncorrelated >/home/rxp163130/scratch/fifth_run/5.2.out" >> uncorrelated
echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 3 -l 5 -p 5 -m uncorrelated >/home/rxp163130/scratch/fifth_run/5.3.out" >> uncorrelated
echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 4 -l 5 -p 10 -m uncorrelated >/home/rxp163130/scratch/fifth_run/5.4.out" >> uncorrelated
echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 5 -l 5 -p 10 -m uncorrelated >/home/rxp163130/scratch/fifth_run/5.5.out" >> uncorrelated

echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 1 -l 62 -p 2 -m uncorrelated >/home/rxp163130/scratch/fifth_run/62.1.out" >> uncorrelated
echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 2 -l 62 -p 3 -m uncorrelated >/home/rxp163130/scratch/fifth_run/62.2.out" >> uncorrelated
echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 3 -l 62 -p 3 -m uncorrelated >/home/rxp163130/scratch/fifth_run/62.3.out" >> uncorrelated
echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 4 -l 62 -p 4 -m uncorrelated >/home/rxp163130/scratch/fifth_run/62.4.out" >> uncorrelated
echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 5 -l 62 -p 4 -m uncorrelated >/home/rxp163130/scratch/fifth_run/62.5.out" >> uncorrelated


echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 1 -l 68 -p 5 -m uncorrelated >/home/rxp163130/scratch/fifth_run/68.1.out" >> uncorrelated
echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 2 -l 68 -p 10 -m uncorrelated >/home/rxp163130/scratch/fifth_run/68.2.out" >> uncorrelated
echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 3 -l 68 -p 10 -m uncorrelated >/home/rxp163130/scratch/fifth_run/68.3.out" >> uncorrelated
echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 4 -l 68 -p 20 -m uncorrelated >/home/rxp163130/scratch/fifth_run/68.4.out" >> uncorrelated
echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 5 -l 68 -p 20 -m uncorrelated >/home/rxp163130/scratch/fifth_run/68.5.out" >> uncorrelated

echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 1 -l 82 -p 1 -m uncorrelated >/home/rxp163130/scratch/fifth_run/82.1.out" >> uncorrelated
echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 2 -l 82 -p 1 -m uncorrelated >/home/rxp163130/scratch/fifth_run/82.2.out" >> uncorrelated
echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 3 -l 82 -p 1 -m uncorrelated >/home/rxp163130/scratch/fifth_run/82.3.out" >> uncorrelated
echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 4 -l 82 -p 1 -m uncorrelated >/home/rxp163130/scratch/fifth_run/82.4.out" >> uncorrelated
echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 5 -l 82 -p 1 -m uncorrelated >/home/rxp163130/scratch/fifth_run/82.5.out" >> uncorrelated

echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 1 -l 99 -p 3 -m uncorrelated >/home/rxp163130/scratch/fifth_run/99.1.out" >> uncorrelated
echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 2 -l 99 -p 4 -m uncorrelated >/home/rxp163130/scratch/fifth_run/99.2.out" >> uncorrelated
echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 3 -l 99 -p 4 -m uncorrelated >/home/rxp163130/scratch/fifth_run/99.3.out" >> uncorrelated
echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 4 -l 99 -p 5 -m uncorrelated >/home/rxp163130/scratch/fifth_run/99.4.out" >> uncorrelated
echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i 5 -l 99 -p 5 -m uncorrelated >/home/rxp163130/scratch/fifth_run/99.5.out" >> uncorrelated

