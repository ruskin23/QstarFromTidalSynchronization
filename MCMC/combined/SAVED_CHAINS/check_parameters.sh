#!/bin/bash

for S in 85 76 96 81 80 36 83 84 94 32 106 123 50 39 56 126 54 70 88 67 95 25 137 1 86 43 73 92 93 79 47 109 44 48 17 8 12 20 31 57 120 28 13
    do
		cat ganymede/MCMC_$S/accepted_parameters_*.txt >ganymede/MCMC_$S/combined_accepted.txt
        cat stampede/MCMC_$S/accepted_parameters_*.txt >stampede/MCMC_$S/combined_accepted.txt
		cat ganymede/MCMC_$S/combined_accepted.txt stampede/MCMC_$S/combined_accepted.txt >combined_accepted.txt
        echo -n "System = $S "
        W="$(wc -l combined_accepted.txt)"
        A="$(cut -d' ' -f1<<<$W)" 
        echo -n "Accepted = ${A} " 
        rm ganymede/MCMC_$S/combined_accepted.txt
        rm stampede/MCMC_$S/combined_accepted.txt
        rm combined_accepted.txt
		
        cat ganymede/MCMC_$S/rejected_parameters_*.txt >ganymede/MCMC_$S/combined_rejected.txt
        cat stampede/MCMC_$S/rejected_parameters_*.txt >stampede/MCMC_$S/combined_rejected.txt
        cat ganymede/MCMC_$S/combined_rejected.txt stampede/MCMC_$S/combined_rejected.txt >combined_rejected.txt
        W="$(wc -l combined_rejected.txt)"
        R="$(cut -d' ' -f1<<<$W)"
        echo -n "Rejected = ${R} " 
        rm ganymede/MCMC_$S/combined_rejected.txt
        rm stampede/MCMC_$S/combined_rejected.txt
        rm combined_rejected.txt
        
        S=$(($A+$R))

        
        r=$(bc <<<"scale=2;$A/$S")
        echo "Acceptance Ratio = $r "
        #printf "\n"
	done | sort -g -k6 

