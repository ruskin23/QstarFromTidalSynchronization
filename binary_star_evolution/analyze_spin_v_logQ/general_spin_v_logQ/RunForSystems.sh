#!/bin/bash

for i in 120 123 125 126 128 129 132 133 134 137 138 142
do
    nohup python3 main.py $i -b 0.0 -n >output/breaks0.0/output_$i.txt 2>&1 &
done
