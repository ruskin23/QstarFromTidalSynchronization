for A in 85 76 96 81 80 36 83 94
do
	for B in 1 2 3 4 5 
	do
		echo "python3 /work/06850/rpatel23/stampede2/QstarFromTidalSynchronization/MCMC/combined/main.py -c -i $B -l $A  -m uncorrelated >"/scratch/06850/rpatel23/output/sixth_run/$A.$B.out
	done
done >NonSyncJob1

for A in 32 106 123 50 39 56 126 54 70 
do
    for B in 1 2 3 4 5
    do
        echo "python3 /work/06850/rpatel23/stampede2/QstarFromTidalSynchronization/MCMC/combined/main.py -c -i $B -l $A  -m uncorrelated >"/scratch/06850/rpatel23/output/sixth_run/$A.$B.out 
    done
done >NonSyncJob2

for A in 88 67 95 25 137 86 43
do
    for B in 1 2 3 4 5
    do
        echo "python3 /work/06850/rpatel23/stampede2/QstarFromTidalSynchronization/MCMC/combined/main.py -c -i $B -l $A  -m uncorrelated >"/scratch/06850/rpatel23/output/sixth_run/$A.$B.out 
    done
done >tmp1
for A in 1 2 3 4 5
do
    echo "python3 /work/06850/rpatel23/stampede2/QstarFromTidalSynchronization/MCMC/combined/main.py -c -i $A -l 1 -m adaptive >"/scratch/06850/rpatel23/output/sixth_run/1.$A.out
done >tmp2
cat tmp1 tmp2 >NonSyncJob3
rm tmp1 tmp2

for A in 73 92 93 79 47 109 44 48 
do
    for B in 1 2 3 4 5
    do
        echo "python3 /work/06850/rpatel23/stampede2/QstarFromTidalSynchronization/MCMC/combined/main.py -c -i $B -l $A  -m uncorrelated >"/scratch/06850/rpatel23/output/sixth_run/$A.$B.out 
    done
done >SyncJob1

for A in 17 8 12 20 31 57 120 28 13
do
    for B in 1 2 3 4 5
    do
        echo "python3 /work/06850/rpatel23/stampede2/QstarFromTidalSynchronization/MCMC/combined/main.py -c -i $B -l $A  -m uncorrelated >"/scratch/06850/rpatel23/output/sixth_run/$A.$B.out 
    done
done >SyncJob2


