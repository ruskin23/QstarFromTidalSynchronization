for A in 1 5 8 12 13 18 20 25 27 28 31 32 33 36 37 39 43 44 47 48 49 50 54 57 59 62 67 68 69 70 73 76 79 80 81 83 84 85 86 88 92 93 94 95 96 99 106 109 111 120 123 125 126 128 129 132 137
do
	for B in 1 2 3 4 5 6 7 8 9 10
	do
		echo "python3 /home/rxp163130/QstarFromTidalSynchronization/MCMC/general/main.py -s -i $B -l $A >output/"'$LAUNCHER_TSK_ID.out' 
	done
done >jobfile
