import os
import sys
import parse_file

#############################
#
# This script processes output from cluster to obtain final set 
# of mutation trajectories that are used for downstream analysis
#
#############################

import os
import sys
import parse_file
from datetime import date

storage_directory='cluster_output_files/'
final_directory=parse_file.data_directory

if len(sys.argv) < 2:
    populations = []
elif sys.argv[1]=='nonmutators':
    populations = parse_file.complete_nonmutator_lines
elif sys.argv[1]=='mutators':
    populations = parse_file.mutator_lines
elif sys.argv[1]=='all':
    populations = parse_file.all_lines
else:
    populations = sys.argv[1:]
    
for population in populations:

    sys.stdout.write("\nProcessing %s...\n" % population)

    # prepare filenames
    merged_timecourse_filename = '%s%s_merged_timecourse.bz2' % (storage_directory, population)
    depth_timecourse_filename = '%s%s_depth_timecourse.bz2' % (storage_directory, population)
    snp_timecourse_filename = '%s%s_snp_timecourse.bz2' % (storage_directory, population)
    indel_timecourse_filename = '%s%s_indel_timecourse.bz2' % (storage_directory, population)
    likelihood_timecourse_filename = '%s%s_likelihood_timecourse.txt' % (storage_directory, population)


    # Filter SNPs and calculate avg depth per sample
    sys.stdout.write('Filtering SNPS and calculating depth...\n')
    return_val = os.system('python filter_snps_and_calculate_depth.py %s %s %s' % (merged_timecourse_filename, depth_timecourse_filename, snp_timecourse_filename))
    if return_val==0:
        sys.stdout.write('Done!\n')
    else:
        sys.stdout.write("Error!\n")
    
    # Call indels 
    sys.stdout.write("Calling indels...\n")
    return_val = os.system('python call_indels.py %s %s' % (merged_timecourse_filename, indel_timecourse_filename))
    if return_val==0:
        sys.stdout.write('Done!\n')
    else:
        sys.stdout.write("Error!\n")
    
    # Annotate pvalues
    sys.stdout.write("Calculating pvalues...\n")
    return_val = os.system('bzcat %s %s %s | ./annotate_pvalues > %s' % (depth_timecourse_filename, snp_timecourse_filename, indel_timecourse_filename, likelihood_timecourse_filename))
    if return_val==0:
        sys.stdout.write('Done!\n')
    else:
        sys.stdout.write("Error!\n")


sys.stdout.write("\n\nTrajectory post-processing output for LTEE metagenomic sequencing project\n")
sys.stdout.write("Date: %s\n\n" % str(date.today()))

os.system('mkdir -p %s' % final_directory)    
    
# Filter and annotate timecourse
sys.stdout.write('Performing final filtering and annotation step...\n')
return_val = os.system('python combined_annotate_timecourse.py %s %s' % (storage_directory, final_directory))
if return_val==0:
    sys.stdout.write('Done!\n')
else:
    sys.stdout.write("Error!\n")
    
# Infer trajectory states in well-mixed HMM
sys.stdout.write('Inferring trajectory states in well-mixed HMM...\n')
return_val = os.system('python calculate_well_mixed_hmm_wrapper.py')
if return_val==0:
    sys.stdout.write('Done!\n')
else:
    sys.stdout.write("Error!\n")

# Infer trajectory states in clade HMM
sys.stdout.write('Inferring trajectory states in clade HMM...\n')
return_val = os.system('python calculate_clade_hmm_wrapper.py')
if return_val==0:
    sys.stdout.write('Done!\n')
else:
    sys.stdout.write("Error!\n")
