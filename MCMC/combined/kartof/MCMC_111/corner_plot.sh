#!/bin/bash


FILE=AcceptedParameters.txt

if test -f "$FILE"; then
    rm $FILE
fi

touch $FILE

for I in 1 2 4 5
do
    python3 FillParameters.py $I 
done

#python3 save_covariance.py $S
python3 corner_plot.py 

#rm temp*.out
