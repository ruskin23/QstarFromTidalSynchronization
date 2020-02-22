#!/bin/bash

for i in 14 18 23 27 33 35 37 49 59
do
    rm ./break0.0/PercentileAges/System_$i/SpinLogQ_10.txt
    #nohup python3 diffages_main.py $i  >output/breaks0.0/PercentileAges/output_$i.txt 2>&1 &
done
