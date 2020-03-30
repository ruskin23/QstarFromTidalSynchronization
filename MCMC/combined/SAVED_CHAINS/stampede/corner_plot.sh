#!/bin/bash

FILE=AcceptedParameters.txt

if test -f "$FILE"; then
    rm $FILE
fi

S=$1
touch $FILE

for I in 3
do
    python3 FillParameters.py $S $I >temp.out
done

python3 corner_plot.py

