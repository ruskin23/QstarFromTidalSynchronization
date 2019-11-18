#!/bin/bash

for i in 13 23 48 82
do
    for j in 0.5 1.0 1.5 2.5 3.0
    do
        nohup python3 main.py $i -b $j >output_files/output_Breaks$j/output_$i.txt 2>&1 &
    done
done
