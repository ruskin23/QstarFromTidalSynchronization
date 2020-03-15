#!/bin/bash


sbatch MCMCNonSync.slurm
sleep 2
sbatch MCMCSync.slurm  
