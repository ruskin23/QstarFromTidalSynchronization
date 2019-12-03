#!/bin/bash

for i in 100 106 109 111 117 119 120 123 125 126 128 129 132 133 134 137 138 142
do
    nohup python3 main.py $i -b 0.5 -n >output/breaks0.5/output_$i.txt 2>&1 &
done
