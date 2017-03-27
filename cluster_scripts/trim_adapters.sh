#!/bin/bash

export population=$1
export sequence_type=$2
mkdir -p data/trimmed_fastq_files
for i in $(python population_parameters.py samples ${population} ${sequence_type}); do
    export sample_name=$i
    export params=$sample_name
    if [ ! -e "data/trimmed_fastq_files/${sample_name}.finished" ]
    then 
    sbatch --job-name=trim-${i} trim_adapters.sbatch
    fi
done
