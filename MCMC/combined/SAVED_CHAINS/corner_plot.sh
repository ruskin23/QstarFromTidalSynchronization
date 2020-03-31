#!/bin/bash


FILE=AcceptedParameters.txt

if test -f "$FILE"; then
    rm $FILE
fi

touch $FILE

S=$1
for I in 1 2 3 4 5
do
    python3 FillParameters.py $S $I 
done

#python3 save_covariance.py $S
python3 corner_plot.py $S

#rm temp*.out
