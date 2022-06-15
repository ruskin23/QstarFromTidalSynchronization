#!/bin/bash

for S in 8 12 13 20 26 28 29 31 33 36 41 43 47 56 63 67 70 76 79 80 82 83 88 90 91 99 109 115 117 120 121 126 128 129 131 139 142
    do
		cat ganymede/MCMC_$S/accepted_parameters_*.txt >combined_accepted.txt
        echo -n "System = $S "
        W="$(wc -l combined_accepted.txt)"
        A="$(cut -d' ' -f1<<<$W)" 
        echo -n "Accepted = ${A} " 
        rm combined_accepted.txt
		
        cat ganymede/MCMC_$S/rejected_parameters_*.txt >combined_rejected.txt
        W="$(wc -l combined_rejected.txt)"
        R="$(cut -d' ' -f1<<<$W)"
        echo -n "Rejected = ${R} " 
        rm combined_rejected.txt
        
        sum=$(($A+$R))

        
        r=$(bc <<<"scale=2;$A/$sum")
        echo "Acceptance Ratio = $r "
        #printf "\n"
	done | sort -g -k6 

