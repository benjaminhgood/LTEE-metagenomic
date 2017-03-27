#!/bin/bash
#SBATCH -p serial_requeue
#SBATCH -t 5:00:00
#SBATCH --mem-per-cpu 16000
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -o trim.out      #File to which standard out will be written
#SBATCH -e trim.err      #File to which standard err will be written

if [[ ${params} == DL_LTE_rn1_* ]]
then
# adaptor trimming
# single ends
java -jar ${TRIMMOMATIC}/trimmomatic-0.32.jar SE -phred33 -threads 1 data/fastq_files/${params}.R0.fastq.gz data/trimmed_fastq_files/${params}.R0.fastq ILLUMINACLIP:${TRIMMOMATIC}/adapters/NexteraPE-PE.fa:2:30:10:2:'false' LEADING:20 TRAILING:20 2> data/trimmed_fastq_files/${params}.trim.unpaired.txt

# paired ends
java -jar ${TRIMMOMATIC}/trimmomatic-0.32.jar PE -phred33 -threads 1 data/fastq_files/${params}.R1.fastq.gz data/fastq_files/${params}.R2.fastq.gz data/trimmed_fastq_files/${params}.R1.fastq data/trimmed_fastq_files/${params}.R1.unpaired.fastq data/trimmed_fastq_files/${params}.R2.fastq data/trimmed_fastq_files/${params}.R2.unpaired.fastq ILLUMINACLIP:${TRIMMOMATIC}/adapters/NexteraPE-PE.fa:2:30:10:2:'false' LEADING:20 TRAILING:20 2> data/trimmed_fastq_files/${params}.trim.paired.txt
fi

if [[ ${params} == ara* ]]
then

java -jar ${TRIMMOMATIC}/trimmomatic-0.32.jar PE -phred33 -threads 1 data/fastq_files/${params}.R1.fastq.gz data/fastq_files/${params}.R2.fastq.gz data/trimmed_fastq_files/${params}.R1.fastq data/trimmed_fastq_files/${params}.R1.unpaired.fastq data/trimmed_fastq_files/${params}.R2.fastq data/trimmed_fastq_files/${params}.R2.unpaired.fastq ILLUMINACLIP:${TRIMMOMATIC}/adapters/NexteraPE-PE.fa:2:30:10:2:'false' LEADING:20 TRAILING:20 2> data/trimmed_fastq_files/${params}.trim.paired.txt

fi

if [[ ${params} == batch* ]]
then

java -jar ${TRIMMOMATIC}/trimmomatic-0.32.jar PE -phred33 -threads 1 data/fastq_files/${params}.R1.fastq.gz data/fastq_files/${params}.R2.fastq.gz data/trimmed_fastq_files/${params}.R1.fastq data/trimmed_fastq_files/${params}.R1.unpaired.fastq data/trimmed_fastq_files/${params}.R2.fastq data/trimmed_fastq_files/${params}.R2.unpaired.fastq ILLUMINACLIP:${TRIMMOMATIC}/adapters/NexteraPE-PE.fa:2:30:10:2:'false' LEADING:20 TRAILING:20 2> data/trimmed_fastq_files/${params}.trim.paired.txt

fi

if [[ ${params} == CL* ]]
then

#lane 1
java -jar ${TRIMMOMATIC}/trimmomatic-0.32.jar PE -phred33 -threads 1 data/fastq_files/${params}.L1.R1.fastq.gz data/fastq_files/${params}.L1.R2.fastq.gz data/trimmed_fastq_files/${params}.L1.R1.fastq data/trimmed_fastq_files/${params}.L1.R1.unpaired.fastq data/trimmed_fastq_files/${params}.L1.R2.fastq data/trimmed_fastq_files/${params}.L1.R2.unpaired.fastq ILLUMINACLIP:${TRIMMOMATIC}/adapters/NexteraPE-PE.fa:2:30:10:2:'false' LEADING:20 TRAILING:20 2> data/trimmed_fastq_files/${params}.L1.trim.paired.txt

#lane 2
java -jar ${TRIMMOMATIC}/trimmomatic-0.32.jar PE -phred33 -threads 1 data/fastq_files/${params}.L2.R1.fastq.gz data/fastq_files/${params}.L2.R2.fastq.gz data/trimmed_fastq_files/${params}.L2.R1.fastq data/trimmed_fastq_files/${params}.L2.R1.unpaired.fastq data/trimmed_fastq_files/${params}.L2.R2.fastq data/trimmed_fastq_files/${params}.L2.R2.unpaired.fastq ILLUMINACLIP:${TRIMMOMATIC}/adapters/NexteraPE-PE.fa:2:30:10:2:'false' LEADING:20 TRAILING:20 2> data/trimmed_fastq_files/${params}.L2.trim.paired.txt

fi

if [[ ${params} == SRR* || ${params} == ERR* ]]
then

if [ ! -e "data/fastq_files/${sample_name}_2.fastq.gz" ]
then
    # single end

java -jar ${TRIMMOMATIC}/trimmomatic-0.32.jar SE -phred33 -threads 1 data/fastq_files/${params}_1.fastq.gz data/trimmed_fastq_files/${params}.R0.fastq ILLUMINACLIP:${TRIMMOMATIC}/adapters/TruSeq3-PE-2.fa:2:30:10:2:'false' LEADING:20 TRAILING:20 2> data/trimmed_fastq_files/${params}.trim.unpaired.txt

else

java -jar ${TRIMMOMATIC}/trimmomatic-0.32.jar PE -phred33 -threads 1 data/fastq_files/${params}_1.fastq.gz data/fastq_files/${params}_2.fastq.gz data/trimmed_fastq_files/${params}.R1.fastq data/trimmed_fastq_files/${params}.R1.unpaired.fastq data/trimmed_fastq_files/${params}.R2.fastq data/trimmed_fastq_files/${params}.R2.unpaired.fastq ILLUMINACLIP:${TRIMMOMATIC}/adapters/TruSeq3-PE-2.fa:2:30:10:2:'false' LEADING:20 TRAILING:20 2> data/trimmed_fastq_files/${params}.trim.paired.txt

fi

fi

touch data/trimmed_fastq_files/${params}.finished
