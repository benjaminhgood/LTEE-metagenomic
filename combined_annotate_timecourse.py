import numpy
import sys
import parse_file
import timecourse_utils

gene_data = parse_file.parse_gene_list()
repeat_data = parse_file.parse_repeat_list()
mask_data = parse_file.parse_mask_list()

position_gene_map, effective_gene_lengths, substitution_specific_synonymous_fraction = parse_file.create_annotation_map(gene_data, repeat_data, mask_data)

input_directory_prefix = sys.argv[1]
output_directory_prefix = sys.argv[2]

FDR = 0.05
pvalue_idx = 6

sys.stdout.write("Filtering trajectories with FDR=%g\n" % FDR)

input_filename_template = input_directory_prefix+"%s_likelihood_timecourse.txt" # change this later
output_filename_template = output_directory_prefix+"%s_annotated_timecourse.txt" 

# first determine pvalue threshold
# across all nonmutator populations
pvalues = []
for population in parse_file.complete_nonmutator_lines:

    file = open(input_filename_template % population,"r") 
    file.readline() # depth line!
    
    for line in file:
        items = line.split(",")
        location = long(items[1])
        pitems = items[pvalue_idx].split()
        pvalue = float(pitems[1])
        autocorrelation = float(pitems[0])
        
        if location==0:
            print line
        
        if parse_file.annotate_gene(location, position_gene_map)=='repeat':
            continue
        
        pvalues.append(pvalue)
    file.close()
    
ntot = len(pvalues)
   
pvalues = numpy.array(pvalues)
   
for p in sorted(pvalues,reverse=True):
    if ntot*p/(pvalues <= p).sum() < FDR:
        break        

# this is now the global threshold pvalue
nonmutator_threshold_pvalue = p
sys.stdout.write("Nonmutator p(FDR) = %g\n" % nonmutator_threshold_pvalue)

# now across all mutator populations
pvalues = []
for population in parse_file.mutator_lines:

    file = open(input_filename_template % population,"r") 
    file.readline() # depth line!
    
    for line in file:
        items = line.split(",")
        location = long(items[1])
        pitems = items[pvalue_idx].split()
        pvalue = float(pitems[1])
        autocorrelation = float(pitems[0])
        
        if location==0:
            print line
        
        if parse_file.annotate_gene(location, position_gene_map)=='repeat':
            continue
        
        pvalues.append(pvalue)
    file.close()
    
ntot = len(pvalues)
   
pvalues = numpy.array(pvalues)
   
for p in sorted(pvalues,reverse=True):
    if ntot*p/(pvalues <= p).sum() < FDR:
        break        

# this is now the global threshold pvalue
mutator_threshold_pvalue = p
sys.stdout.write("Mutator p(FDR) = %g\n" % mutator_threshold_pvalue)

total_num_passed = 0

# now go through all the populations and annotate stuff based on this
for population in parse_file.all_lines:

    if population in parse_file.complete_nonmutator_lines:
        threshold_pvalue = nonmutator_threshold_pvalue
    else:
        threshold_pvalue = mutator_threshold_pvalue

    # total number of trajectories processed
    num_total = 0
    # total passed
    num_passed = 0
    
    annotated_mutations = []
    header = None
    printed_header = False
    loaded_avg_depths = False

    file = open(input_filename_template % population,"r") 
    for line in file:
        
        num_total += 1
    
        # parse timecourse data
        items = line.split(",")
        chromosome = items[0].strip()
        location = long(items[1])
        allele = items[2].strip()
        total_times = numpy.array([float(subitem) for subitem in items[3].split()])*1000
        total_alts = numpy.array([float(subitem) for subitem in items[4].split()])
        total_depths = numpy.array([float(subitem) for subitem in items[5].split()])
        test_statistic = float(items[pvalue_idx].split()[0])
        pvalue = float(items[pvalue_idx].split()[1])
        deletion_idx, fold_reduction, deletion_pvalue = tuple([float(subitem) for subitem in items[pvalue_idx+1].split()])
        deletion_idx = long(deletion_idx)
        
        duplication_idx, fold_increase, duplication_pvalue = tuple([float(subitem) for subitem in items[pvalue_idx+2].split()])
        duplication_idx = long(duplication_idx)
    
        times = total_times[total_times<1000000]
        alts = total_alts[total_times<1000000]
        depths = total_depths[total_times<1000000]
        
        clone_times = total_times[total_times>1000000]-1000000
        clone_alts = total_alts[total_times>1000000]
        clone_depths = total_depths[total_times>1000000]
        
        # load average depths if not leaded yet
        if not loaded_avg_depths:
            pop_avg_depths = depths
            clone_avg_depths = clone_depths
            loaded_avg_depths = True
        
        # create header if not created yet
        if header == None:    
            print_strings = ['Position', 'Gene', 'Allele', 'Annotation', 'Test statistic', 'P-value', 'Deletion index', 'Fold reduction', 'Deletion P-value', 'Duplication index', 'Fold increase', 'Duplication pvalue', 'Passed?']
            for t in zip(total_times):
                print_strings.append('AC:%d' % t)
                print_strings.append('DP:%d' % t)
            
            header = ", ".join(print_strings)
        
        # annotate mutation
        gene_name, var_type = parse_file.annotate_variant(location, allele, gene_data, position_gene_map)
        
        # if pvalue is lower than threshold and not a weird insertion gene
        passed_str = 'FAIL'
        
        if (pvalue <= threshold_pvalue) and gene_name!='repeat' and (location!=parse_file.ancestral_araA_mutation_location) and (location!=parse_file.ancestral_recD_mutation_location) and (alts[0]*1.0/(depths[0]+(depths[0]==0)) < 0.2):
            # would otherwise pass
            # check if clone fail
            
            # determine whether clone data suggests that the mutation 
            # is a duplication. do not pass these
        
            # first estimate frequencies at good timepoints    
            good_idxs, filtered_alts, filtered_depths = timecourse_utils.mask_timepoints(times, alts, depths, var_type, deletion_idx, fold_reduction, deletion_pvalue)
            freqs = timecourse_utils.estimate_frequencies(filtered_alts, filtered_depths)
            masked_times = times[good_idxs]
            masked_freqs = freqs[good_idxs]
        
            masked_depth_ratios = depths[good_idxs]/pop_avg_depths[good_idxs]
        
            interpolation_function = timecourse_utils.create_interpolation_function(masked_times, masked_freqs, tmax=100000)
    
            masked_clone_times, masked_clone_freqs = timecourse_utils.estimate_clone_frequencies(clone_times, clone_alts, clone_depths)
  
            if len(masked_clone_times) == 0:
                # clone fail.. there are no points with actual reads.. this is an issue
                clone_fail = True
            
            else:
            
                # we have some clone timepoints to work with
                clone_fail = False
                
                # try to detect duplication without clone information
                if (masked_freqs>0.25).sum() > 10: # sticks around for a while
                    if (masked_depth_ratios[masked_freqs>0.25].mean()/(masked_depth_ratios[:10].mean()) > 1.5): # evidence of duplication
                        if (((masked_freqs>0.25)*(masked_freqs<0.75)).sum()*1.0/(masked_freqs>0.25).sum() > 0.9): # never really fixes
                            if (masked_clone_freqs<0.9).all(): # just to make sure we aren't picking a trajectory that has fixed in a clone
                                clone_fail = True
            
                # now try to detect deletion with clone information
    
                # frequency in population at time the clone is called
                pop_freqs = interpolation_function(masked_clone_times)
    
                
                # number of timepoints where mutation is at intermediate frequency in a clone
                intermediate_clone_idxs = ((masked_clone_freqs>0.10)*(masked_clone_freqs<0.7))
                num_intermediate_clones = intermediate_clone_idxs.sum()
                if intermediate_clone_idxs.sum() >= 4:
                    if ((pop_freqs>0.2)*(pop_freqs<0.7)*intermediate_clone_idxs).sum()*1.0/intermediate_clone_idxs.sum() >= 0.5:
                        clone_fail = True
        
                if ((masked_clone_freqs<0.6)*(pop_freqs<0.6)).all() and ((masked_clone_freqs>0.1)*(pop_freqs>0.1)).any():
                    clone_fail = True
        
                # see if there is evidence for a duplication
                # calculate clone depth changes specifically where coverage is ok
                clone_freqs = clone_alts*1.0/(clone_depths+(clone_depths==0))
                clone_depth_fold_changes = timecourse_utils.estimate_depth_fold_changes(clone_avg_depths, clone_depths)
                clone_idxs = numpy.array([t in masked_clone_times for t in clone_times])
                clone_depth_fold_changes = clone_depth_fold_changes[clone_idxs] 
                num_duplicated_intermediates = (clone_depth_fold_changes[(masked_clone_freqs>0.4)*(masked_clone_freqs<0.8)]>0.4).sum()*1.0
    
                if ( (num_duplicated_intermediates > 0) and (masked_clone_freqs<0.9).all() and (masked_freqs<0.9).all() ):
                    clone_fail = True
    
                if masked_freqs[0]>0.1:
                    clone_fail=True
                
            if not clone_fail:
                passed_str = 'PASS'
                num_passed += 1
        
        # print to CSV file
        print_strings = [str(location), gene_name, allele, var_type, str(test_statistic), str(pvalue), str(deletion_idx), str(fold_reduction), str(deletion_pvalue), str(duplication_idx), str(fold_increase), str(duplication_pvalue), passed_str]
        for alt,depth in zip(total_alts, total_depths):
            print_strings.append(str(alt))
            print_strings.append(str(depth))
        
        annotated_mutations.append((location, ", ".join(print_strings)))
        
    file.close()

    # now write them out in sorted order
    output_file = open(output_filename_template % population,"w")
    output_file.write(header)
    output_file.write("\n")
    for location, mutation_str in sorted(annotated_mutations,key=lambda x: x[0]):
        output_file.write(mutation_str)
        output_file.write("\n")
    output_file.close()   

    sys.stdout.write("%s: %d passed / %d\n" % (population, num_passed, num_total))

    total_num_passed += num_passed

# done iterating through populations
sys.stdout.write("%d total mutations!\n" % total_num_passed)