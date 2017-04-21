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

from numpy.random import shuffle

import figure_utils
import timecourse_utils
import stats_utils

mpl.rcParams['font.size'] = 6
mpl.rcParams['lines.linewidth'] = 1.0
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'

####
#
# Set up figures
#
####

fig = plt.figure(figsize=(6.2, 1.7))

grid = gridspec.GridSpec(1, 2, width_ratios=[1,1], wspace=0.1)

time_axis = plt.Subplot(fig, grid[0])
fig.add_subplot(time_axis)

time_axis.spines['top'].set_visible(False)
time_axis.spines['right'].set_visible(False)
time_axis.get_xaxis().tick_bottom()
time_axis.get_yaxis().tick_left()

time_axis.set_ylabel('Fraction mutations $\geq t$')
time_axis.set_xlabel('Appearance time, $t$')
time_axis.set_xticks(figure_utils.time_xticks)
time_axis.set_xticklabels(figure_utils.time_xticklabels)
time_axis.set_xlim([0,60000])

missense_time_axis = plt.Subplot(fig, grid[1])
fig.add_subplot(missense_time_axis)


missense_time_axis.spines['top'].set_visible(False)
missense_time_axis.spines['right'].set_visible(False)
missense_time_axis.get_xaxis().tick_bottom()
missense_time_axis.get_yaxis().tick_left()

missense_time_axis.set_xlabel('Appearance time, $t$')
missense_time_axis.set_xticks(figure_utils.time_xticks)
missense_time_axis.set_xticklabels(figure_utils.time_xticklabels)
missense_time_axis.set_xlim([0,60000])
missense_time_axis.set_yticklabels([])

####
#
# Do calculation
#
####

excluded_types = set(['sv'])

# map of variant types 
appearance_times = {}

pooled_appearance_times = []
pooled_var_types = []
restricted_appearance_times = []
restricted_var_types = []

populations = parse_file.complete_nonmutator_lines
        
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
 
        num_processed_mutations+=1 
            
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
        
        pooled_appearance_times.append(t)
        pooled_var_types.append(var_type)
        
        if var_type not in excluded_types:
            restricted_appearance_times.append(t)
            restricted_var_types.append(var_type)
        
    sys.stderr.write("added %d mutations\n" % num_processed_mutations)
    

# sort the lists (but preserve association between time and variant type
pooled_appearance_times, pooled_var_types = (list(x) for x in zip(*sorted(zip(pooled_appearance_times, pooled_var_types), key=lambda pair: pair[0])))
restricted_appearance_times, restricted_var_types = (list(x) for x in zip(*sorted(zip(restricted_appearance_times, restricted_var_types), key=lambda pair: pair[0])))

pooled_appearance_times = numpy.array(pooled_appearance_times)
restricted_appearance_times = numpy.array(restricted_appearance_times)

# calculate appearance time distribution for each variant type
observed_appearance_times = {}
observed_restricted_appearance_times = {}

for t,var_type in zip(pooled_appearance_times, pooled_var_types):
    if var_type not in observed_appearance_times:
        observed_appearance_times[var_type] = []
    observed_appearance_times[var_type].append(t)

for var_type in observed_appearance_times.keys():
    observed_appearance_times[var_type] = numpy.array(observed_appearance_times[var_type])
    observed_appearance_times[var_type].sort()  

for t,var_type in zip(restricted_appearance_times, restricted_var_types):
    if var_type not in observed_restricted_appearance_times:
        observed_restricted_appearance_times[var_type] = []
    observed_restricted_appearance_times[var_type].append(t)
    
for var_type in observed_restricted_appearance_times.keys():
    observed_restricted_appearance_times[var_type] = numpy.array(observed_restricted_appearance_times[var_type])
    observed_restricted_appearance_times[var_type].sort()  

####
#
# Now do bootstrap resampling to gauge significance
#
####

num_bootstraps = 10000
bootstrapped_kss = {var_type:[] for var_type in observed_appearance_times.keys()}
bootstrapped_restricted_kss = {var_type:[] for var_type in observed_restricted_appearance_times.keys()}

sys.stderr.write('Resampling %d bootstraps...\n' % num_bootstraps)
for bootstrap_idx in xrange(0,num_bootstraps):

    # permute var type labels
    shuffle(pooled_var_types)
    shuffle(restricted_var_types) 
    
    # recalculate distributions for each type
    bootstrapped_appearance_times = {}
    bootstrapped_restricted_appearance_times = {}

    for t,var_type in zip(pooled_appearance_times, pooled_var_types):
        if var_type not in bootstrapped_appearance_times:
            bootstrapped_appearance_times[var_type] = []
        bootstrapped_appearance_times[var_type].append(t)

    for var_type in bootstrapped_appearance_times.keys():
        bootstrapped_appearance_times[var_type].sort()  

    for t,var_type in zip(restricted_appearance_times, restricted_var_types):
        if var_type not in bootstrapped_restricted_appearance_times:
            bootstrapped_restricted_appearance_times[var_type] = []
        bootstrapped_restricted_appearance_times[var_type].append(t)
    
    for var_type in bootstrapped_restricted_appearance_times.keys():
        bootstrapped_restricted_appearance_times[var_type].sort()  

    # recalculate ks distances
    for var_type in bootstrapped_appearance_times.keys():
        D = stats_utils.calculate_ks_distance(bootstrapped_appearance_times[var_type], pooled_appearance_times)
        bootstrapped_kss[var_type].append(D)
        
    for var_type in bootstrapped_restricted_appearance_times.keys():
        D = stats_utils.calculate_ks_distance(bootstrapped_restricted_appearance_times[var_type], restricted_appearance_times)
        bootstrapped_restricted_kss[var_type].append(D)

# calculate pvalues
sys.stdout.write('Pvalues for temporal nonuniformity (n=%d bootstraps):\n' % num_bootstraps)
for var_type in observed_appearance_times.keys():

    bootstrapped_appearance_times[var_type] = numpy.array(bootstrapped_appearance_times[var_type])

    D = stats_utils.calculate_ks_distance(observed_appearance_times[var_type], pooled_appearance_times)
    
    pvalue = ((bootstrapped_kss[var_type]>=D).sum()+1.0)/(len(bootstrapped_kss[var_type])+1.0)
        
    sys.stdout.write('%s: %g\n' % (var_type, pvalue))

sys.stdout.write('Excluding svs:\n')
for var_type in observed_restricted_appearance_times.keys():

    bootstrapped_restricted_appearance_times[var_type] = numpy.array(bootstrapped_restricted_appearance_times[var_type])

    D = stats_utils.calculate_ks_distance(observed_restricted_appearance_times[var_type], restricted_appearance_times)
    
    pvalue = ((bootstrapped_restricted_kss[var_type]>=D).sum()+1.0)/(len(bootstrapped_restricted_kss[var_type])+1.0)
        
    sys.stdout.write('%s: %g\n' % (var_type, pvalue))
    
######
#
# Now do plotting
#
######            

all_ts, all_survivals = stats_utils.calculate_unnormalized_survival_from_vector(pooled_appearance_times, min_x=-1000,max_x=100000)

time_axis.step(all_ts, all_survivals/all_survivals[0], color='k', label='All')
missense_time_axis.step(all_ts, all_survivals/all_survivals[0], color='k', label='All')

restricted_ts, restricted_survivals = stats_utils.calculate_unnormalized_survival_from_vector(restricted_appearance_times, min_x=-1000,max_x=100000)
missense_time_axis.step(all_ts, restricted_survivals/restricted_survivals[0], color='k', label='All (excluding sv)',alpha=0.5)



for var_type in parse_file.var_types:
 
    color = figure_utils.get_var_type_color(var_type)
    vartype_ts, vartype_survivals = stats_utils.calculate_unnormalized_survival_from_vector(observed_appearance_times[var_type], min_x=-1000, max_x=100000)
    time_axis.step(vartype_ts, vartype_survivals/vartype_survivals[0], color=color, alpha=0.7, label=var_type)

    if var_type == 'missense':
        missense_time_axis.step(vartype_ts, vartype_survivals/vartype_survivals[0], color=color, alpha=0.7, label=var_type)
    

time_axis.legend(loc='upper right', frameon=False)
missense_time_axis.legend(loc='upper right', frameon=False)

sys.stderr.write("Saving figure...\t")
fig.savefig(parse_file.figure_directory+'supplemental_vartype_times.pdf',bbox_inches='tight')
sys.stderr.write("Done!\n")