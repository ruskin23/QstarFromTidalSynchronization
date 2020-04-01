#!/bin/bash



nohup python3 diffages_main.py $1 -p $2 >output/breaks0.0/PercentileAges/$1.$2.txt 2>&1 &
