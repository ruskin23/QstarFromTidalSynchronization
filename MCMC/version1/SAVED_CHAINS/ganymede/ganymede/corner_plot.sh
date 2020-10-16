#!/bin/bash

FILE=AcceptedParameters.txt

if test -f "$FILE"; then
    rm $FILE
fi

S=$1
touch $FILE

for I in 1 2 3 4 5
do
    python3 FillParameters.py $S $I >temp.out
done

python3 corner_plot.py

