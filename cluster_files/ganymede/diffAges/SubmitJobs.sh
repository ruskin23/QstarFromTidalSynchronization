#!/bin/bash

sbatch Ages_1.slurm  
sleep 2
sbatch Ages_2.slurm
sleep 2
sbatch Ages_3.slurm
