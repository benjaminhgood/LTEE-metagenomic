#!/bin/bash

export population=$1
export sequence_type=$2
export gd_files="" 

for i in $(python population_parameters.py samples ${population} ${sequence_type}); do
    export sample_name=$i

    if [ ! -e "data/breseq_output/${sample_name}/output/evidence/evidence.gd" ]
    then
    echo "Should not get here!"
    rm data/breseq_output/${sample_name}/output/filtered_output.gd
    fi

    echo ${sample_name}
    cat data/breseq_output/${sample_name}/output/evidence/evidence.gd | grep 'JC\|#' > data/breseq_output/${sample_name}/output/filtered_output.gd 

    export gd_files="${gd_files} data/breseq_output/${sample_name}/output/filtered_output.gd"
    
done

gdtools UNION -o data/breseq_output/${population}_merged_breseq_output.gd -e ${gd_files}
