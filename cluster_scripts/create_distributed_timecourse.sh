#!/bin/bash
export population=$1
export sequence_type=$2
export bam_files=""
export sample_names=""
#export reference_file=data/lenski_reference/REL606.fasta
export reference_file=data/rebreseq_output/${population}_DL_LTE_rn1_96_GTAGAGGA-CTAAGCCT/data/reference.fasta

for i in $(python population_parameters.py samples ${population} ${sequence_type}); do
    
    export sample_name=$i
    export sample_names="${sample_names} ${sample_name}"
    export bam_files="${bam_files} data/rebreseq_output/${population}_${sample_name}/data/reference.bam"
    #export bam_files="${bam_files} data/bam_files/${sample_name}_sorted_dedupped_readgroups.bam"
done
#echo ${bam_files}

#rm data/timecourse_files/*
#for i in {50000..50000..50000}
for i in {50000..4700000..50000}
do
    export start_position=`expr $i - 50000`
    export end_position=$i
    #bash create_timecourse.sbatch
    sbatch --job-name=timecourse_${population}_${start_position} create_timecourse.sbatch
done
