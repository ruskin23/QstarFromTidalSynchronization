for A in 1 8 12 13 25 32 36 39 43 47 50 54 56 67 70 76 80 81 83 84 85 86 88 94 95 96 106 109 123 126 137
do
	for B in 1 2 3 4 
	do
		echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i $B -l $A >/home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/output/ganymede/"'$LAUNCHER_TSK_ID.out' 
	done
done >jobfile_2
