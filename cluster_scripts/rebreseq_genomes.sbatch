#!/bin/bash
#SBATCH -p general
#SBATCH -t 108:00:00
#SBATCH --mem 48000
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -o rebreseq.out      #File to which standard out will be written
#SBATCH -e rebreseq.err      #File to which standard err will be written
# old mem was 16000
# old mem was 11750
# highgest mem was 32000

mkdir -p data/rebreseq_output/${population}_${sample_name}
rm -rf data/rebreseq_output/${population}_${sample_name}
mkdir -p data/rebreseq_output/${population}_${sample_name}

export FASTQ_FILES=""
for f in $(ls data/trimmed_fastq_files/${sample_name}*.fastq)
do
   export FASTQ_FILES="${FASTQ_FILES} ${f}"
done

breseq -p -j 1 --user-evidence-gd data/breseq_output/${population}_merged_breseq_output.gd -o data/rebreseq_output/${population}_${sample_name}/ -r data/lenski_reference/REL606.6.gbk ${FASTQ_FILES} > data/rebreseq_output/${population}_${sample_name}.out 2> data/rebreseq_output/${population}_${sample_name}.err
