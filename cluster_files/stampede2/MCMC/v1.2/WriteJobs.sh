for A in 12 47 106 109 8 28 25 17 83  
do
	for B in 1 2 3 4 5 
	do
		echo "python3 /work/06850/rpatel23/stampede2/QstarFromTidalSynchronization/MCMC/version1/main.py -c -i $B -l $A  -m uncorrelated >"/scratch/06850/rpatel23/output/Run44/$A.$B.out
	done
done >uncorrelated1

for A in  31 79 120 44 94 13 48 43 93 
do
    for B in 1 2 3 4 5
    do
        echo "python3 /work/06850/rpatel23/stampede2/QstarFromTidalSynchronization/MCMC/version1/main.py -c -i $B -l $A  -m uncorrelated >"/scratch/06850/rpatel23/output/Run44/$A.$B.out 
    done
done >uncorrelated2

for A in 57 36 88 70
do
    for B in 1 2 3 4 5
    do
        echo "python3 /work/06850/rpatel23/stampede2/QstarFromTidalSynchronization/MCMC/version1/main.py -c -i $B -l $A  -m uncorrelated >"/scratch/06850/rpatel23/output/Run44/$A.$B.out
    done
done >mixed

for A in 50
do
    for B in 1 2 3 4 5
    do
        echo "python3 /work/06850/rpatel23/stampede2/QstarFromTidalSynchronization/MCMC/version1/main.py -c -i $B -l $A  -m adaptive >"/scratch/06850/rpatel23/output/Run44/$A.$B.out 
    done
done >>mixed



for A in 39 137 54 80 126 76 1 32 67
do
    for B in 1 2 3 4 5
    do
        echo "python3 /work/06850/rpatel23/stampede2/QstarFromTidalSynchronization/MCMC/version1/main.py -c -i $B -l $A  -m adaptive >"/scratch/06850/rpatel23/output/Run44/$A.$B.out 
    done
done >adaptive1

for A in 81 95 96 85 56 73 86 92 20 
do
    for B in 1 2 3 4 5
    do
        echo "python3 /work/06850/rpatel23/stampede2/QstarFromTidalSynchronization/MCMC/version1/main.py -c -i $B -l $A  -m adaptive >"/scratch/06850/rpatel23/output/Run44/$A.$B.out 
    done
done >adaptive2

#done >adaptive3

