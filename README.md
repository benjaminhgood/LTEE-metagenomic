LICENSE: The code located here is provided as-is, with no warranty, etc. (It's under GPL v2.) 

REQUIREMENTS: Python (v2.7) with the NumPy (v1.10.4) and SciPy (v0.17.1) libraries, and a C++ compiler that supports the C++11 standard (e.g., GNU's g++). 

INSTALLATION: (1) Compile the "main.cpp" source file in the "trajectory_pvalue_cpp_code" directory, name the program "annotate_pvalues", and move it to the base directory. 

USAGE: (1) To generate figures, run "python process_cluster_output.py", followed by "python generate_all_figures.py".  

AUXILLIARY PROGRAMS: For completeness, we have also included the scripts used to process the raw FASTQ files to obtain an initial set of candidate mutations. These scripts are designed to run on a computing cluster with the slurm job queuing system, and require trimmomatic (v0.32), samtools, and python/numpy/scipy to be installed in the local environment. We have also included a modified version of the breseq pipeline (Deatherage and Barrick, Methods Mol Biol, 2014) named breseq-lite. This program is copyright its respective authors and subject to other licenses as detailed in the open-source files distributed with breseq.

JOURNAL REFERENCE: BH Good*, MJ McDonald*, JE Barrick, RE Lenski, and MM Desai. The Dynamics of Molecular Evolution Over 60,000 Generations. Nature, 551:45-50 (2017). 

