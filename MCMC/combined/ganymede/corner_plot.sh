#!/bin/bash

FILE=AcceptedParameters.txt

if test -f "$FILE"; then
    rm $FILE
fi

S=$1
touch $FILE

for I in 1 2 3 4 5
do
	echo "Instance = $I"
    python3 FillParameters.py $S $I 
done

cat MCMC_$1/rejected_parameters_* >combined_rejected.txt
cat MCMC_$1/accepted_parameters_* >combined_accepted.txt

wc -l combined_rejected.txt
wc -l combined_accepted.txt
rm combined_rejected.txt
rm combined_accepted.txt

#python3 corner_plot.py

