#!/bin/bash

export population=$1
export sequence_type=$2
mkdir -p data/rebreseq_output

if [ -e "data/breseq_output/${population}_merged_breseq_output.gd" ]
then

for i in $(python population_parameters.py samples ${population} ${sequence_type}); do
    export sample_name=$i
    export params=$sample_name
    
    if [ ! -e "data/rebreseq_output/${population}_${sample_name}/output/evidence/evidence.gd" ]
    then

        sbatch --job-name=rebreseq-${population}_${sample_name} rebreseq_genomes.sbatch
    fi
done

fi
