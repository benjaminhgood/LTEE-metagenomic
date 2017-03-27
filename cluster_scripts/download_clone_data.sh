#!/bin/bash

export index=0
for i in $(python population_parameters.py clone_accessions); do
    index=$((index+1))
    accessions[index]=$i
done

export accessions

#sbatch --array=1-${index} download_clone_data.sbatch

for idx in 63; do
# 127 128 138 140 141 146 150 155 156 159 163
    export SLURM_ARRAY_TASK_ID=$idx
    sbatch download_clone_data.sbatch
done
