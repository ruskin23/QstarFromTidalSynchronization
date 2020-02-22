for B in 1 2 3 4 5
do
	for A in 85 76 96 81 80 36 83 84 94 32 106 123 50 47 39 56 126 54 109 70 8 12 88 67 95 25 137 1 86 43 13
	do
		echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/main.py -s -i $B -l $A >/home/rxp163130/QstarFromTidalSynchronization/MCMC/combined/output/ganymede/"'$LAUNCHER_TSK_ID.out' 
	done
done >jobfile
