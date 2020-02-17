for A in 85 76 96 81 80 36 83 84 94 
do
	for B in 1 2 3 4 5 6 7 8 9 10 
	do
		echo "python3 /work/06850/rpatel23/stampede2/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i $B -l $A >"/scratch/06850/rpatel23/output/'$LAUNCHER_TSK_ID.out' 
	done
done >JobFile
