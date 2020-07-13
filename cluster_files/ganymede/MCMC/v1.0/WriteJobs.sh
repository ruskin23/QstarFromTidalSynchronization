

for B in 1 2 3 4 5
do
    echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -c -i $B -l 1 -m adaptive >/home/rxp163130/scratch/fifth_run/1.$B.out"
	for A in 85 76 96 81 80 36 83 94 32 106 123 50 39 56 126 54 70 88 67 95 25 137 86 43 
	do
		echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -c -i $B -l $A -m uncorrelated >/home/rxp163130/scratch/fifth_run/$A.$B.out"
	done
done >NonSyncJob



for B in 1 2 3 4 5
do
    for A in 73 92 93 79 47 109 44 48 17 8 12 20 31 57 120 28 13 
    do
        echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -c -i $B -l $A > -m uncorrelated /home/rxp163130/scratch/fifth_run/$A.$B.out"
    done
done >SyncJob

