for A in 17 25 109 120 83 57 36 88
do
	for B in 1 2 3 4 5 
	do
		echo "python3 /work/06850/rpatel23/stampede2/QstarFromTidalSynchronization/MCMC/combined/main.py -c -i $B -l $A  -m uncorrelated >"/scratch/06850/rpatel23/output/fourteenth_run/$A.$B.out
	done
done >uncorrelated1

for A in 44 47 79 94 13 48 43 
do
    for B in 1 2 3 4 5
    do
        echo "python3 /work/06850/rpatel23/stampede2/QstarFromTidalSynchronization/MCMC/combined/main.py -c -i $B -l $A  -m uncorrelated >"/scratch/06850/rpatel23/output/fourteenth_run/$A.$B.out 
    done
done >uncorrelated2

for A in 93 70 8 
do
    for B in 1 2 3 4 5
    do
        echo "python3 /work/06850/rpatel23/stampede2/QstarFromTidalSynchronization/MCMC/combined/main.py -c -i $B -l $A  -m uncorrelated >"/scratch/06850/rpatel23/output/fourteenth_run/$A.$B.out 
    done
done >common

for A in 39 137 54 80 126 76 1 32 67
do
    for B in 1 2 3 4 5
    do
        echo "python3 /work/06850/rpatel23/stampede2/QstarFromTidalSynchronization/MCMC/combined/main.py -c -i $B -l $A  -m adaptive >"/scratch/06850/rpatel23/output/fourteenth_run/$A.$B.out 
    done
done >adaptive1

for A in 81 95 96 50 85 56 73 86 92 
do
    for B in 1 2 3 4 5
    do
        echo "python3 /work/06850/rpatel23/stampede2/QstarFromTidalSynchronization/MCMC/combined/main.py -c -i $B -l $A  -m adaptive >"/scratch/06850/rpatel23/output/fourteenth_run/$A.$B.out 
    done
done >adaptive2

for A in 20
do
    for B in 1 2 3 4 5
    do
        echo "python3 /work/06850/rpatel23/stampede2/QstarFromTidalSynchronization/MCMC/combined/main.py -c -i $B -l $A  -m adaptive >"/scratch/06850/rpatel23/output/fourteenth_run/$A.$B.out  >> common        
    done
done
#done >adaptive3

