import pylab
import numpy
import sys
from math import log10
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import matplotlib.patheffects as pe

from numpy.random import binomial
from scipy.interpolate import interp1d
from scipy.special import btdtr
import bz2
import parse_file

import matplotlib
import matplotlib.pyplot as plt

import timecourse_utils
import figure_utils

mutator_color=figure_utils.mutator_group_color
nonmutator_color=figure_utils.nonmutator_group_color

min_fixation_coverage = 10

mpl.rcParams['font.size'] = 5.0
mpl.rcParams['lines.linewidth'] = 1.0
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'

###########################################
# Set up figures
###########################################

##############################################################################
#
# Supplemental fig (comparing minority and majority fixation trajectories)
#
##############################################################################

minority_fig = plt.figure(figsize=(2.2, 1.5))
minority_grid = gridspec.GridSpec(1, 1)
minority_axis = plt.Subplot(minority_fig, minority_grid[0])

minority_fig.add_subplot(minority_axis)

minority_axis.set_xlabel('Generation, $t$')
minority_axis.set_ylabel('Num fixed mutations')

minority_axis.spines['top'].set_visible(False)
minority_axis.spines['right'].set_visible(False)
minority_axis.get_xaxis().tick_bottom()
minority_axis.get_yaxis().tick_left()

minority_axis.set_xlim([-1000,61000])
minority_axis.set_ylim([0,100])
minority_axis.set_xticks(figure_utils.time_xticks)
minority_axis.set_xticklabels(figure_utils.time_xticklabels)

majority_fixed_states = set([parse_file.clade_hmm_states['FB'], parse_file.clade_hmm_states['FM'], parse_file.clade_hmm_states['PB*']])
minority_fixed_states = set([parse_file.clade_hmm_states['FB'], parse_file.clade_hmm_states['Fm'], parse_file.clade_hmm_states['PB*']])

######
#
# Now do the calculations
#
######

populations = ['p5','m6']

theory_times = numpy.arange(0,121)*500.0
for population in populations:

    majority_Ms = numpy.zeros_like(theory_times)
    minority_Ms = numpy.zeros_like(theory_times)

    # load mutations
    mutations, depth_tuple = parse_file.parse_annotated_timecourse(population)
    population_avg_depth_times, population_avg_depths, clone_avg_depth_times, clone_avg_depths = depth_tuple
    dummy_times,fmajors,fminors,haplotype_trajectories = parse_file.parse_haplotype_timecourse(population)
    state_times, state_trajectories = parse_file.parse_well_mixed_state_timecourse(population)
    
    for mutation_idx in xrange(0,len(mutations)):

        location, gene_name, allele, var_type, test_statistic, pvalue, cutoff_idx, depth_fold_change, depth_change_pvalue, times, alts, depths, clone_times, clone_alts, clone_depths = mutations[mutation_idx] 

        Ls = haplotype_trajectories[mutation_idx]
        state_Ls = state_trajectories[mutation_idx]
        
        good_idxs, filtered_alts, filtered_depths = timecourse_utils.mask_timepoints(times, alts, depths, var_type, cutoff_idx, depth_fold_change, depth_change_pvalue)
    
        freqs = timecourse_utils.estimate_frequencies(filtered_alts, filtered_depths)
        
        masked_times = times[good_idxs]
        masked_freqs = freqs[good_idxs]
        masked_depths = depths[good_idxs]
        masked_state_Ls = state_Ls[good_idxs]
        masked_Ls = Ls[good_idxs]    
        
        t0,tf,dt = timecourse_utils.calculate_appearance_fixation_time_from_clade_hmm(masked_times, masked_freqs, masked_Ls)
        
        if masked_Ls[-1] in majority_fixed_states:
            majority_Ms += (theory_times>tf)   
        if masked_Ls[-1] in minority_fixed_states:
            minority_Ms += (theory_times>tf)    
        
        
    color = parse_file.get_line_color(population)
    
    
    line, = minority_axis.plot(theory_times, minority_Ms, ':', linewidth=0.5, color=color)
    line.set_dashes((1,0.5))
    minority_axis.plot(theory_times, majority_Ms, '-', linewidth=0.5, color=color)
    
minority_axis.plot([-1],[-1],'k-',linewidth=0.5,label='Major')
line, = minority_axis.plot([-1],[-1],'k:',linewidth=0.5,label='Minor')
line.set_dashes((1,0.5))
    
minority_axis.legend(loc='upper left',frameon=False)
    
    

############
#
# Save figure
#
############

sys.stderr.write("Saving figure...\t")
minority_fig.savefig(parse_file.figure_directory+'supplemental_majority_minority_comparison.pdf', bbox_inches='tight')
sys.stderr.write("Done!\n")