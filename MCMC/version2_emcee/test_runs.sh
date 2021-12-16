#!/bin/bash

for S in {1..20}
do nohup python3 mcmc_sampler.py >output_files/output_$S.txt 2>&1 &
done
