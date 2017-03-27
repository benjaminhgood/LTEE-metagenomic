#!/bin/bash

export population=$1
export sequence_type=$2

for i in $(python population_parameters.py samples ${population} ${sequence_type}); do
    export sample_name=$i
    export params=$sample_name
  
    python lines_exceed.py 10000 data/fastq_files/${sample_name}*.fastq.gz
 
done
