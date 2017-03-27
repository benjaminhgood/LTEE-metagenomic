#!/bin/bash
#SBATCH -p serial_requeue
#SBATCH -t 20:00:00
#SBATCH --mem-per-cpu 8000
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -o download.out      #File to which standard out will be written
#SBATCH -e download.err      #File to which standard err will be written

export accession=$(python population_parameters.py clone_accession_index ${SLURM_ARRAY_TASK_ID})

cd data/fastq_files
fastq-dump -I --split-files $accession
