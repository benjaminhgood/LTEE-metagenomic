######
#
# This script plots Fig S6 
# (distribution of appearance times of different variant types)
#
#####

import sys
import numpy
import parse_file
from scipy.special import gammaln as loggamma
from math import log, exp

import matplotlib as mpl
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patheffects as pe
from mpl_toolkits.axes_grid.inset_locator import inset_axes

from numpy.random import shuffle, choice, normal

import figure_utils
import timecourse_utils
import stats_utils

mpl.rcParams['font.size'] = 6
mpl.rcParams['lines.linewidth'] = 1.0
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'

bar_width=0.35

metapopulations = ['nonmutator','mutator']
metapopulation_labels = {'nonmutator': 'Nonmutators', 'mutator':'Mutators'}
metapopulation_populations = {'nonmutator': parse_file.complete_nonmutator_lines, 'mutator': parse_file.mutator_lines}
metapopulation_offsets = {'nonmutator': -bar_width, 'mutator': 0}
metapopulation_colors = {'nonmutator': figure_utils.nonmutator_group_color, 'mutator': figure_utils.mutator_group_color}

####
#
# Set up figure
#
####

fig = plt.figure(figsize=(4, 1.7))

grid = gridspec.GridSpec(1, 2)

dnds_axis = plt.Subplot(fig, grid[0])
fig.add_subplot(dnds_axis)

dnds_axis.spines['top'].set_visible(False)
dnds_axis.spines['right'].set_visible(False)
dnds_axis.get_xaxis().tick_bottom()
dnds_axis.get_yaxis().tick_left()

dnds_axis.set_ylabel('dN/dS')
dnds_axis.set_ylim([0,3.5])
dnds_axis.set_yticks([0,1,2,3,4])
dnds_axis.set_xlim([-1,2])
dnds_axis.set_xticks([0,1])
dnds_axis.set_xticklabels(['All','Fixed'],rotation='vertical')
dnds_axis.plot([-4,4],[1,1],'-',linewidth=0.25,color='0.7')

dnds_axis.set_title('Naive')

targeted_dnds_axis = plt.Subplot(fig, grid[1])
fig.add_subplot(targeted_dnds_axis)

targeted_dnds_axis.spines['top'].set_visible(False)
targeted_dnds_axis.spines['right'].set_visible(False)
targeted_dnds_axis.get_xaxis().tick_bottom()
targeted_dnds_axis.get_yaxis().tick_left()

targeted_dnds_axis.set_ylim([0,3.5])
targeted_dnds_axis.set_yticks([0,1,2,3,4])
targeted_dnds_axis.set_xlim([-1,2])
targeted_dnds_axis.set_xticks([0,1])
targeted_dnds_axis.set_xticklabels(['All','Fixed'],rotation='vertical')
targeted_dnds_axis.plot([-4,4],[1,1],'-',linewidth=0.25,color='0.7')

targeted_dnds_axis.set_title('Substitution specific')
####
#
# Do calculation
#
####

nonsynonymous_types = set(['missense','nonsense'])
synonymous_types = set(['synonymous'])

Lsyn, Lnon, substitution_specific_synonymous_fraction = parse_file.calculate_synonymous_nonsynonymous_target_sizes()


for metapopulation in metapopulations:

    populations = metapopulation_populations[metapopulation]
    
    non_appeared = {population: 0.0 for population in populations}
    non_fixed = {population: 0.0 for population in populations}
    
    syn_appeared = {population: 0.0 for population in populations}
    syn_fixed = {population: 0.0 for population in populations}
    
    targeted_Lsyn = {population: 0.0 for population in populations}
    targeted_Lnon = {population: 0.0 for population in populations}
    targeted_fixed_Lsyn = {population: 0.0 for population in populations}
    targeted_fixed_Lnon = {population: 0.0 for population in populations}
    
    
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
        
            good_idxs, filtered_alts, filtered_depths = timecourse_utils.mask_timepoints(times, alts, depths, var_type, cutoff_idx, depth_fold_change, depth_change_pvalue)
    
            freqs = timecourse_utils.estimate_frequencies(filtered_alts, filtered_depths)
        
            masked_times = times[good_idxs]
            masked_freqs = freqs[good_idxs]
            masked_state_Ls = state_Ls[good_idxs]
            masked_Ls = Ls[good_idxs]
            
            t = timecourse_utils.calculate_appearance_time(masked_times, masked_freqs, masked_state_Ls, masked_Ls)
            #fixed_weight = timecourse_utils.calculate_fixed_weight(masked_state_Ls[-1], masked_freqs[-1])
            fixed_weight = timecourse_utils.calculate_clade_fixed_weight(masked_Ls[-1], masked_freqs[-1])
        
            
            if var_type in nonsynonymous_types or var_type in synonymous_types: 
                targeted_Lnon[population] += (1-substitution_specific_synonymous_fraction[allele])  
                targeted_fixed_Lnon[population] += fixed_weight*(1-substitution_specific_synonymous_fraction[allele])
                targeted_Lsyn[population] += substitution_specific_synonymous_fraction[allele]
                targeted_fixed_Lsyn[population] += fixed_weight*substitution_specific_synonymous_fraction[allele]
            
            if var_type in nonsynonymous_types:    
                non_appeared[population]+=1
                non_fixed[population]+=fixed_weight
                num_processed_mutations+=1     
                
            elif var_type in synonymous_types:
                syn_appeared[population]+=1
                syn_fixed[population]+=fixed_weight
                num_processed_mutations+=1 
                
                
            
            else:
                pass
            
        sys.stderr.write("added %d mutations\n" % num_processed_mutations)
    
    total_non_appeared = sum([non_appeared[population] for population in populations])
    total_non_fixed = sum([non_fixed[population] for population in populations])
    total_syn_appeared = sum([syn_appeared[population] for population in populations])
    total_syn_fixed = sum([syn_fixed[population] for population in populations])
    
    
    dnds_appeared = total_non_appeared/total_syn_appeared*Lsyn/Lnon 
    dnds_fixed = total_non_fixed/total_syn_fixed*Lsyn/Lnon
    
    total_targeted_Lnon = sum(targeted_Lnon.values())
    total_targeted_Lsyn = sum(targeted_Lsyn.values())
    targeted_dnds_appeared = total_non_appeared/total_syn_appeared*total_targeted_Lsyn/total_targeted_Lnon 
    total_targeted_Lnon = sum(targeted_fixed_Lnon.values())
    total_targeted_Lsyn = sum(targeted_fixed_Lsyn.values())
    targeted_dnds_fixed = total_non_fixed/total_syn_fixed*total_targeted_Lsyn/total_targeted_Lnon
    
    population_dnds_appeared = numpy.array([non_appeared[population]/(syn_appeared[population]+(syn_appeared[population]==0))*Lsyn/Lnon for population in populations])
    
    population_dnds_fixed = numpy.array([non_fixed[population]/(syn_fixed[population]+(syn_fixed[population]==0))*Lsyn/Lnon for population in populations])
    
    targeted_population_dnds_appeared = numpy.array([non_appeared[population]/syn_appeared[population]*targeted_Lsyn[population]/targeted_Lnon[population] for population in populations])
    
    targeted_population_dnds_fixed = numpy.array([non_fixed[population]/(syn_fixed[population]+(syn_fixed[population]==0))*targeted_fixed_Lsyn[population]/targeted_fixed_Lnon[population] for population in populations])
    
    
    for population in populations:
        sys.stdout.write("%s: %d non, %d syn, dN/dS = %g, (dN/dS)* = %g\n" % (population, non_appeared[population], syn_appeared[population], non_appeared[population]/syn_appeared[population]*Lsyn/Lnon, non_appeared[population]/syn_appeared[population]*targeted_Lsyn[population]/targeted_Lnon[population]))
        sys.stdout.write("%s fixed: %d non, %d syn\n" % (population, non_fixed[population], syn_fixed[population]))
        
    sys.stdout.write("Total: %d non, %d syn, dN/dS = %g, (dN/dS)* = %g\n" % (total_non_appeared, total_syn_appeared, dnds_appeared, targeted_dnds_appeared))
        
    sys.stderr.write("Bootstrap resampling dNdS...\t")

    bootstrapped_dnds_appeared = []
    bootstrapped_dnds_fixed = []
    num_bootstraps = 10000
    for bootstrap_idx in xrange(0,num_bootstraps):
        bootstrapped_populations = choice(populations, len(populations))
        
        bootstrapped_non_appeared = sum([non_appeared[population] for population in bootstrapped_populations])
        bootstrapped_non_fixed = sum([non_fixed[population] for population in bootstrapped_populations])
        bootstrapped_syn_appeared = sum([syn_appeared[population] for population in bootstrapped_populations])
        bootstrapped_syn_fixed = sum([syn_fixed[population] for population in bootstrapped_populations])
    
        bootstrapped_dnds_appeared.append( bootstrapped_non_appeared/bootstrapped_syn_appeared*Lsyn/Lnon )
        bootstrapped_dnds_fixed.append( bootstrapped_non_fixed/bootstrapped_syn_fixed*Lsyn/Lnon )
    
    bootstrapped_dnds_appeared = numpy.array( bootstrapped_dnds_appeared )
    bootstrapped_dnds_fixed = numpy.array( bootstrapped_dnds_fixed )
    sys.stderr.write("Done!\n")
    
    sys.stderr.write("Calculating confidence intervals...\t")
    bootstrapped_dnds_appeared.sort()
    bootstrapped_dnds_fixed.sort()
    
    upper_dnds_appeared = stats_utils.calculate_percentile_from_sorted_vector(bootstrapped_dnds_appeared, 0.975)
    lower_dnds_appeared = stats_utils.calculate_percentile_from_sorted_vector(bootstrapped_dnds_appeared, 0.025)
    
    upper_dnds_fixed = stats_utils.calculate_percentile_from_sorted_vector(bootstrapped_dnds_fixed, 0.975)
    lower_dnds_fixed = stats_utils.calculate_percentile_from_sorted_vector(bootstrapped_dnds_fixed, 0.025)
    sys.stderr.write("Done!\n")    
 
    ######
    #
    # Now do plotting
    #
    ######  
    
    dnds_axis.bar([metapopulation_offsets[metapopulation], 1+metapopulation_offsets[metapopulation]],[dnds_appeared, dnds_fixed],width=bar_width, bottom=-0.05, linewidth=0, facecolor=metapopulation_colors[metapopulation],alpha=0.5)
    
    dnds_axis.plot( numpy.ones_like(population_dnds_appeared)*(metapopulation_offsets[metapopulation]+bar_width/2)+normal(0,bar_width/2*0.3*numpy.ones_like(population_dnds_appeared)),population_dnds_appeared,'.',color=metapopulation_colors[metapopulation],alpha=0.5,label=metapopulation_labels[metapopulation])
    dnds_axis.plot( numpy.ones_like(population_dnds_fixed)*(metapopulation_offsets[metapopulation]+bar_width/2+1)+normal(0,bar_width/2*0.3*numpy.ones_like(population_dnds_fixed)),population_dnds_fixed,'.',color=metapopulation_colors[metapopulation],alpha=0.5)
    
    targeted_dnds_axis.bar([metapopulation_offsets[metapopulation], 1+metapopulation_offsets[metapopulation]],[targeted_dnds_appeared, targeted_dnds_fixed],width=bar_width, bottom=-0.05, linewidth=0, facecolor=metapopulation_colors[metapopulation],alpha=0.5)
    
    targeted_dnds_axis.plot( numpy.ones_like(targeted_population_dnds_appeared)*(metapopulation_offsets[metapopulation]+bar_width/2)+normal(0,bar_width/2*0.3*numpy.ones_like(targeted_population_dnds_appeared)),targeted_population_dnds_appeared,'.',color=metapopulation_colors[metapopulation],alpha=0.5)
    targeted_dnds_axis.plot( numpy.ones_like(targeted_population_dnds_fixed)*(metapopulation_offsets[metapopulation]+bar_width/2+1)+normal(0,bar_width/2*0.3*numpy.ones_like(targeted_population_dnds_fixed)),targeted_population_dnds_fixed,'.',color=metapopulation_colors[metapopulation],alpha=0.5)
   
dnds_axis.legend(loc='upper right', frameon=False) 
sys.stderr.write("Saving figure...\t")          
fig.savefig(parse_file.figure_directory+'extended_data_fig3.pdf',bbox_inches='tight')
sys.stderr.write("Done!\n")