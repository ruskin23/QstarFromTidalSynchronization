#!/bin/bash

for i in 14 17 23 35 56 63 65 78 82 100 117 119 133 134 138 142
do
    nohup python3 main.py $i -b 0.0 -a >output/output_$i.txt 2>&1 &
done
