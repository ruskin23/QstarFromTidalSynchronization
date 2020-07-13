#!/bin/bash

FILE=extracted.txt

if test -f "$FILE"; then
    rm $FILE
fi

touch $FILE


for S in 85 76 96 81 80 36 83 84 94 32 106 123 50 39 56 126 54 70 88 67 95 25 137 1 86 43 73 92 93 79 47 109 44 48 17 8 12 20 57 120 28 13
    do python3 extract.py $S
done

