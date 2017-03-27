#!/bin/bash
export population=$1
export sequence_type=$2
export gd_files=""
export sample_names=""
export reference_file=data/lenski_reference/REL606.fasta

mkdir -p data/breseq_timecourse_files/

for i in $(python population_parameters.py samples ${population} ${sequence_type}); do
    export sample_name=$i
    export sample_names="${sample_names} ${sample_name}"
    export gd_files="${gd_files} data/rebreseq_output/${population}_${sample_name}/output/evidence/evidence.gd"
done

# write everything to a merged file
cat data/timecourse_files/${population}_*timecourse.txt | python create_breseq_timecourse.py ${reference_file} ${population} ${gd_files} | bzip2 -c > ${population}_merged_timecourse.bz2
