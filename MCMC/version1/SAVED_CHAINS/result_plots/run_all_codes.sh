#!/bin/bash

nohup python3 detail_fitting.py >temp_1.out 2>&1 &
sleep 3600
nohup python3 p_values_class.py p_E >temp_2.out 2>&1 &
sleep 3600
nohup python3 p_value_class.py E_p >temp_3.out 2>&1 &
sleep 3600
nohup python3 plotting.py >temp_4.txt 2>&1 &

