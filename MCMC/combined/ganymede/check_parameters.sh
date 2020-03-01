#!/bin/bash

for S in 85 76 96 81 80 36 83 84 94 32 106 123 50 39 56 126 54 70 88 67 95 25 137 1 86 43 73 92 93 79 47 109 44 48 17 8 12 20 31 57 120 28 13
	do
		cat MCMC_$S/accepted_parameters_*.txt >MCMC_$S/combined_accepted.txt
		wc -l MCMC_$S/combined_accepted.txt
		rm MCMC_$S/combined_accepted.txt
        cat MCMC_$S/rejected_parameters_*.txt >MCMC_$S/combined_rejected.txt
		wc -l MCMC_$S/combined_rejected.txt
        rm MCMC_$S/combined_rejected.txt
		printf "\n"
	done
