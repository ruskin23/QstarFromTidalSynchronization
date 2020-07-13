#!/bin/bash

for S in 111 128 132 18 27 33 49 5 62 68 99
    do
        echo -n "System = $S "
        W="$(wc -l samples/MassAgeFehSamples_$S.txt)"
        A="$(cut -d' ' -f1<<<$W)" 
        echo -n "Accepted = ${A} " 
	
        W="$(wc -l rejected_parameters/rejected_parameters_$S.txt)"
        R="$(cut -d' ' -f1<<<$W)"
        echo -n "Rejected = ${R} " 
        
        diff=$(($A-$R))
        
        echo "Actual Accepted = $diff"
        
        r=$(bc <<<"scale=2;$diff/$A")
        echo "Acceptance Ratio = $r "
        #printf "\n"
	done 

