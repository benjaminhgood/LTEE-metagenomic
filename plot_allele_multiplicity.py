####
#
# This script plots Fig. S7, which looks at patterns of allele multiplicity 
#
###

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
import figure_utils

mpl.rcParams['font.size'] = 6
mpl.rcParams['lines.linewidth'] = 1.0
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'

####
#
# Set up figure
#
####

fig = plt.figure(figsize=(2.75, 3.25))

grid = gridspec.GridSpec(2, 1, height_ratios=[1,1], hspace=0.15)

nonmutator_axis = plt.Subplot(fig, grid[0])
fig.add_subplot(nonmutator_axis)

nonmutator_axis.spines['top'].set_visible(False)
nonmutator_axis.spines['right'].set_visible(False)
nonmutator_axis.get_xaxis().tick_bottom()
nonmutator_axis.get_yaxis().tick_left()

mutator_axis = plt.Subplot(fig, grid[1])
fig.add_subplot(mutator_axis)

mutator_axis.set_ylabel('Total mutations $\geq m$',y=1.1)
mutator_axis.set_xlabel('Nucleotide multiplicity, $m$')

mutator_axis.spines['top'].set_visible(False)
mutator_axis.spines['right'].set_visible(False)
mutator_axis.get_xaxis().tick_bottom()
mutator_axis.get_yaxis().tick_left()

####
#
# Do calculation
#
####

excluded_types = set(['sv','indel'])

reference_sequence = parse_file.parse_reference_genome()
gene_data = parse_file.parse_gene_list()
repeat_data = parse_file.parse_repeat_list()
mask_data = parse_file.parse_mask_list()
    
position_gene_map, effective_gene_lengths, substitution_specific_synonymous_fraction = parse_file.create_annotation_map(gene_data, repeat_data, mask_data)

#Ltot = 4.4e06
Ltot = len(reference_sequence)-effective_gene_lengths['masked']
sys.stderr.write("Ltot = %d\n" % Ltot)    

for population_group in ['nonmutators','mutators']:

    if population_group == 'nonmutators':
        populations = parse_file.complete_nonmutator_lines
        color = figure_utils.nonmutator_group_color
        label = figure_utils.nonmutator_group_label
        axis = nonmutator_axis
        yrange = [1,1e03]
        
    elif population_group == 'mutators':
        populations = parse_file.mutator_lines
        color = figure_utils.mutator_group_color
        label = figure_utils.mutator_group_label
        axis = mutator_axis
        yrange = [1,1e05]
        
    else:
        continue
        
    position_population_map = {}
    
    multiplicities = []
        
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
        
            if var_type not in excluded_types:
                
                num_processed_mutations+=1
                
                if position not in position_population_map:
                    position_population_map[position] = set()
                    
                position_population_map[position].add(population)
         
        sys.stderr.write("added %d mutations\n" % num_processed_mutations)
    
    sys.stdout.write("Sites mutated >= 6 times:\n")            
    for position in position_population_map.keys():
        m = len(position_population_map[position])
        if m>=6:
            sys.stderr.write("%d\n" % position)
        multiplicities.append( len(position_population_map[position]) )
        
    multiplicities.sort()
    multiplicities = numpy.array(multiplicities)*1.0
    ntot = multiplicities.sum()
    
    ks = numpy.arange(1,8)
    observed_survival = numpy.array([multiplicities[(multiplicities>=k)].sum() for k in ks])
    observed_survival = numpy.clip(observed_survival, 1e-01, 1000000)
    
    # Null distribution (SI) 
    poisson_pmf = numpy.array([exp((k-1)*log(ntot/Ltot)-loggamma(k)-ntot/Ltot) for k in ks])
    poisson_survival = numpy.array([poisson_pmf[(k-1):].sum() for k in ks])
    null_survival = ntot*poisson_survival
    
    
    axis.step(ks, null_survival, where='mid',color='0.7',linewidth=0.5)
    axis.step(ks, observed_survival, where='mid', color=color)
    axis.semilogy([0],[1])
    axis.set_xlim([0.4,6.6])
    axis.set_ylim(yrange)
    sys.stdout.write("%s: %d mutations in >=2 hit sites (%g of total)\n" % (label, observed_survival[1], observed_survival[1]/observed_survival[0]))
    
nonmutator_axis.set_xticklabels([])
sys.stderr.write("Saving figure...\t")
fig.savefig(parse_file.figure_directory+'extended_data_fig4.pdf',bbox_inches='tight') 
sys.stderr.write("Done!\n")    