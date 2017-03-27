Brief instructions:

Required programs: trimmomatic, breseq-lite, and samtools must be installed on the local path.

Data: the raw FASTQ files must be stored in a directory ./data/fastq_files/ (note: data can be a symlink to somewhere else, e.g. temporary storage). You must have write access to the ./data directory. The metadata for the fastq files are stored in the population_samples.csv (and clone_samples.csv) files. You must modify these as necessary for your own experiment. 

Parallelization: the scripts are designed to be run with the SLURM job manager. If you wish to run them locally, simply change "sbatch" commands in the *.sh scripts to "bash" commands, and everything should still work. 

Instructions:

(1) Trim adapters. For each population_id in population_samples.csv, run "bash trim_adapters.sh population_id".  

(2) Initial breseq round. For each population_id in population_samples.csv,
run "bash breseq_genomes.sh population_id". 

(3) Compile merged list of candidate junctions per population. For each population_id in population_samples.csv, run "bash merge_candidate_junctions.sh population_id".

(4) Second round of breseq. For each population_id in population_samples.csv,
run "bash rebreseq_genomes.sh population_id".

(5) Create SNP timecourses. For each population_id in population_samples.csv,
run "bash create_distributed_timecourse.sh population_id".

(6) Create junction timecourses, and final merged timecourse file per
population. For each population_id in population_samples.csv, run "bash
create_merged_timecourse.sh population_id"

End result: popluation_id_merged_timecourse.bz2, which must be processed with
additional scripts.  
