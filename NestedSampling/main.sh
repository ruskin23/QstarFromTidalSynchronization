#!/bin/bash


nohup python3 main.py -l $1 >output_test_$1.txt 2>&1 &
