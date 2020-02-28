for A in 85 76 96 81 80 36 83 84 94
do
	for B in 1 2 3 4 5 
	do
		echo "python3 /work/06850/rpatel23/stampede2/QstarFromTidalSynchronization/MCMC/combined/main.py -c -i $B -l $A >"/scratch/06850/rpatel23/output/second_run/$A.$B.out
	done
done >NonSyncJob1

for A in 32 106 123 50 39 56 126 54 70 
do
    for B in 1 2 3 4 5
    do
        echo "python3 /work/06850/rpatel23/stampede2/QstarFromTidalSynchronization/MCMC/combined/main.py -c -i $B -l $A >"/scratch/06850/rpatel23/output/second_run/$A.$B.out 
    done
done >NonSyncJob2

for A in 88 67 95 25 137 1 86 43
do
    for B in 1 2 3 4 5
    do
        echo "python3 /work/06850/rpatel23/stampede2/QstarFromTidalSynchronization/MCMC/combined/main.py -c -i $B -l $A >"/scratch/06850/rpatel23/output/second_run/$A.$B.out 
    done
done >NonSyncJob3


for A in 73 92 93 79 47 109 44 48 
do
    for B in 1 2 3 4 5
    do
        echo "python3 /work/06850/rpatel23/stampede2/QstarFromTidalSynchronization/MCMC/combined/main.py -c -i $B -l $A >"/scratch/06850/rpatel23/output/second_run/$A.$B.out 
    done
done >SyncJob1

for A in 17 8 12 20 31 57 120 28 13
do
    for B in 1 2 3 4 5
    do
        echo "python3 /work/06850/rpatel23/stampede2/QstarFromTidalSynchronization/MCMC/combined/main.py -c -i $B -l $A >"/scratch/06850/rpatel23/output/second_run/$A.$B.out 
    done
done >SyncJob2


