import sys
import numpy
import parse_file
import timecourse_utils
from numpy.random import shuffle
import stats_utils

metapopulations = ['nonmutator','mutator']
metapopulation_populations = {'nonmutator': parse_file.complete_nonmutator_lines, 'mutator': parse_file.mutator_lines}


def calculate_test_statistic(pfixes):
    return (pfixes-pfixes[0])[1:].min()
    #return numpy.diff(pfixes).min()

for metapopulation in metapopulations:

    populations = metapopulation_populations[metapopulation]

    if metapopulation == 'nonmutator':       
        desired_var_types = ['synonymous','missense', 'nonsense','noncoding','indel','sv']
    else:
        desired_var_types = ['synonymous','missense','nonsense']
    
    var_types = []
    fixed_weights = []

    for population in populations:
    
        sys.stderr.write("Processing %s...\t" % parse_file.get_pretty_name(population))

        # calculate mutation trajectories
        # load mutations
        mutations, depth_tuple = parse_file.parse_annotated_timecourse(population)
    
        population_avg_depth_times, population_avg_depths, clone_avg_depth_times, clone_avg_depths = depth_tuple
    
        dummy_times,fmajors,fminors,haplotype_trajectories = parse_file.parse_haplotype_timecourse(population)
        state_times, state_trajectories = parse_file.parse_well_mixed_state_timecourse(population)
    
        num_processed_mutations = 0
        for mutation_idx in xrange(0,len(mutations)):
 
            position, gene_name, allele, var_type, test_statistic, pvalue, cutoff_idx, depth_fold_change, depth_change_pvalue, times, alts, depths, clone_times, clone_alts, clone_depths = mutations[mutation_idx] 
            Ls = haplotype_trajectories[mutation_idx]
            state_Ls = state_trajectories[mutation_idx]
        
            if var_type not in desired_var_types:
                continue
                
            num_processed_mutations+=1 
            
            good_idxs, filtered_alts, filtered_depths = timecourse_utils.mask_timepoints(times, alts, depths, var_type, cutoff_idx, depth_fold_change, depth_change_pvalue)
    
            freqs = timecourse_utils.estimate_frequencies(filtered_alts, filtered_depths)
        
            masked_times = times[good_idxs]
            masked_freqs = freqs[good_idxs]
            masked_state_Ls = state_Ls[good_idxs]
            masked_Ls = Ls[good_idxs]
            
            t = timecourse_utils.calculate_appearance_time(masked_times, masked_freqs, masked_state_Ls, masked_Ls)
            #fixed_weight = timecourse_utils.calculate_fixed_weight(masked_state_Ls[-1], masked_freqs[-1])
            fixed_weight = timecourse_utils.calculate_clade_fixed_weight(masked_Ls[-1], masked_freqs[-1])
        
            
            var_types.append(var_type)
            fixed_weights.append(fixed_weight)
            
        sys.stderr.write("added %d mutations!\n" % num_processed_mutations)

    # Calculate pfixes
    num_total = {var_type: 0 for var_type in desired_var_types}
    num_fixed = {var_type: 0 for var_type in desired_var_types}
    for var_type, fixed_weight in zip(var_types, fixed_weights):
        num_total[var_type]+=1.0
        num_fixed[var_type]+=fixed_weight
    
    pfixes = numpy.array([num_fixed[var_type]/num_total[var_type] for var_type in desired_var_types])

    observed_delta = calculate_test_statistic(pfixes)
    
    print pfixes, observed_delta
    
    sys.stderr.write("Bootstrapping fixation probabilities...\t")
    if metapopulation=='nonmutator':
        num_bootstraps = 10000
    else:
        num_bootstraps = 1000
    bootstrapped_deltas = []
    for bootstrap_idx in xrange(0,num_bootstraps):
        
        shuffle(fixed_weights)
        
        bootstrapped_num_total = {var_type: 0 for var_type in desired_var_types}
        bootstrapped_num_fixed = {var_type: 0 for var_type in desired_var_types}
        for var_type, fixed_weight in zip(var_types, fixed_weights):
            bootstrapped_num_total[var_type]+=1.0
            bootstrapped_num_fixed[var_type]+=fixed_weight
    
        bootstrapped_pfixes = numpy.array([bootstrapped_num_fixed[var_type]/bootstrapped_num_total[var_type] for var_type in desired_var_types])
        
        bootstrapped_delta = calculate_test_statistic(bootstrapped_pfixes)
        
        if bootstrapped_delta > observed_delta:
            pass
            #print bootstrapped_pfixes
            #print bootstrapped_delta_pfixes
        
            
        bootstrapped_deltas.append( bootstrapped_delta )
    
    bootstrapped_deltas = numpy.array(bootstrapped_deltas)
    sys.stderr.write("Done!\n")
    
    pvalue = stats_utils.calculate_empirical_pvalue(observed_delta, bootstrapped_deltas)
    
    sys.stdout.write("Survival probability sum test for %ss, p=%g (n=%d bootstraps)\n" % (metapopulation, pvalue, num_bootstraps))
    