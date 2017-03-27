############
#
# Calculates convergence matrices for different mutation identity classes
#
# These are of the form: 
#
############

import pylab
import numpy
import sys
from math import log10
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
from numpy.random import binomial
from scipy.interpolate import interp1d
import bz2
import parse_file
from scipy.special import gammaln
from math import exp,log

import matplotlib
import matplotlib.pyplot as plt

import timecourse_utils

total_times = numpy.hstack([numpy.arange(0,60000.0/500+1)*500.0])

###################
#
# Load gene information
#
###################

sys.stderr.write("Calculating gene sizes...\t")
gene_size_map = parse_file.create_gene_size_map()
sys.stderr.write("Done!\n")

###################
#
# Load operon information
#
###################

sys.stderr.write("Loading operon information...\t")
operon_data = parse_file.parse_operon_list()
operon_size_map = operon_data[2]
sys.stderr.write("Done!\n")

# Make four different types of convergence matrices
# Gene level, including svs
# Gene level, no svs
# Operon level, ...
matrix_types = [('gene', True), ('gene', False), ('operon',True), ('operon',False)]

populations = parse_file.all_lines

for level, include_svs in matrix_types:

    sys.stderr.write('Calculating convergence matrix at %s level (svs=%r)...\n' % (level, include_svs))

    excluded_types = set(['synonymous'])
    if not include_svs:
        excluded_types.add('sv')
    
    if level=='gene':
        identifier_size_map = gene_size_map
        
        def get_identifier_name(gene_name):
            return gene_name
        
    elif level=='operon':
        identifier_size_map = operon_size_map
        
        def get_identifier_name(gene_name):
            return parse_file.annotate_operon(gene_name, operon_data)
    else:
        sys.stderr.write("Should never get here!\n")
        
    convergence_matrix = {}
    for identifier_name in sorted(identifier_size_map.keys()):
        #print identifier_name
        
        length = identifier_size_map[identifier_name]
        #length = max([length,200])
        
        convergence_matrix[identifier_name] = {'length': length, 'mutations': {population: [] for population in populations}}
            
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
        
            if gene_name=='intergenic':
                continue
                
            if var_type in excluded_types:
                continue
        
            identifier = get_identifier_name(gene_name)
            
            if identifier==None:
                sys.stderr.write("No identifier for %s!\n" % gene_name)
                continue
            
            num_processed_mutations += 1    
            
            good_idxs, filtered_alts, filtered_depths = timecourse_utils.mask_timepoints(times, alts, depths, var_type, cutoff_idx, depth_fold_change, depth_change_pvalue)
    
            freqs = timecourse_utils.estimate_frequencies(filtered_alts, filtered_depths)
        
            masked_times = times[good_idxs]
            masked_freqs = freqs[good_idxs]
            masked_state_Ls = state_Ls[good_idxs]
            masked_Ls = Ls[good_idxs]
            
            t = timecourse_utils.calculate_appearance_time(masked_times, masked_freqs, masked_state_Ls, masked_Ls)
        
            convergence_matrix[identifier]['mutations'][population].append((t, masked_state_Ls[-1], masked_Ls[-1], masked_freqs[-1]))

        sys.stderr.write("processed %d mutations!\n" % num_processed_mutations)

    # Print it out
    if include_svs:
        output_filename = parse_file.data_directory+("%s_convergence_matrix.txt" % level)
    else:
        output_filename = parse_file.data_directory+("%s_convergence_matrix_nosvs.txt" % level)
        
    convergence_matrix_file = open(output_filename,"w")
    
    
    # Header
    convergence_matrix_file.write(", ".join(["Identifier"]+["Size"]+[population for population in populations]))
    
    for identifier in sorted(convergence_matrix.keys()):

        length = convergence_matrix[identifier]['length']
        mutations = convergence_matrix[identifier]['mutations']

        convergence_matrix_file.write("\n")
        convergence_matrix_file.write(", ".join([identifier, "%0.1f" % length]+[";".join(["%d:%d:%d:%g" % (t,L,Lclade,f) for t,L,Lclade,f in mutations[population]]) for population in populations]))
    
    convergence_matrix_file.close()