#!/bin/bash
#oldSBATCH -p serial_requeue
#SBATCH -p general
#oldSBATCH -t 47:00:00
#oldSBATCH --mem 11750
#oldSBATCH -t 72:00:00
#oldSBATCH --mem 16000
#SBATCH -t 108:00:00
#SBATCH --mem 20000
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -o breseq.out      #File to which standard out will be written
#SBATCH -e breseq.err      #File to which standard err will be written

export FASTQ_FILES=""
for f in $(ls data/trimmed_fastq_files/${sample_name}*.fastq)
do
   export FASTQ_FILES="${FASTQ_FILES} ${f}"
done

breseq -p --brief-html-output --polymorphism-reject-indel-homopolymer-length 0 --polymorphism-reject-surrounding-homopolymer-length 0 --polymorphism-score-cutoff 2 --polymorphism-minimum-coverage-each-strand 0 -o data/breseq_output/${sample_name} -r data/lenski_reference/REL606.6.gbk ${FASTQ_FILES} > data/breseq_output/${sample_name}.out 2> data/breseq_output/${sample_name}.err

