#!/bin/bash

#2 5 14 18 23 27 33 35 37 49 59 62 63 65 68 69 77 78 82 99 100 111 117 119 125 128 129 132 133 134 138 142

for  S in 18 
do
	for  P in 1
		do python3 diffages_main.py $S -p $P >output/breaks0.0/PercentileAges/$S.$P.txt 2>&1 &
	
    done
done 

