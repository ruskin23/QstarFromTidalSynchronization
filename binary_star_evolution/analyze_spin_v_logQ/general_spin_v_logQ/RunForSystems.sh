#!/bin/bash

for i in 62 69 59 129 5 111 99 18 27 33 49 132 128 37 125 68
do
    nohup nohup python3 diffages_main.py $i -n >output/breaks0.0/PercentileAges/output_$i.txt 2>&1 &
done
