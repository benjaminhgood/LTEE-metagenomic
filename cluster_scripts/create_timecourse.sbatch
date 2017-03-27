#!/bin/bash
#SBATCH -p serial_requeue
#SBATCH -t 8:00:00
#SBATCH --mem 3700
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -o create_full_timecourse.out      #File to which standard out will be written
#SBATCH -e create_full_timecourse.err      #File to which standard err will be written


#old one for pre-breseq pipeline
#samtools mpileup -q10 -f ${reference_file} -r chr1:${start_position}-${end_position} ${bam_files} > data/timecourse_files/${population}_${start_position}.pileup
samtools mpileup -q10 -f ${reference_file} -r REL606:${start_position}-${end_position} ${bam_files} > data/timecourse_files/${population}_${start_position}.pileup

cat data/timecourse_files/${population}_${start_position}.pileup | python create_timecourse.py ${start_position} ${end_position} ${population} ${sample_names} > data/timecourse_files/${population}_${start_position}_timecourse.txt

# remove temporary file
# UNCOMMENT TO SAVE SPACE
rm data/timecourse_files/${population}_${start_position}.pileup

#bzcat realigned_timecourse.bz2 | python split_timecourse.py ${start_position} ${end_position} > data/timecourse_files/${start_position}_timecourse.txt 
