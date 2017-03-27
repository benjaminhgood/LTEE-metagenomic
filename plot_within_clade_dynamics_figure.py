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

mpl.rcParams['font.size'] = 5.0
mpl.rcParams['lines.linewidth'] = 1.0
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'

###########################################
# Set up figures
###########################################

##############################################################################
#
# First set up Fig. 4 (Within-clade dynamics)
#
##############################################################################

fig = plt.figure(figsize=(5, 2.75))


# make three panels panels
outer_grid  = gridspec.GridSpec(1, 2, width_ratios=[1,1], wspace=0.30)

inner_grid_1 = gridspec.GridSpecFromSubplotSpec(2, 1, height_ratios=[1,1],
                subplot_spec=outer_grid[1], hspace=0.15) #, hspace=0.08)

inner_grid_2 = gridspec.GridSpecFromSubplotSpec(2, 1, height_ratios=[1,1],
                subplot_spec=outer_grid[0], hspace=0.15) #, hspace=0.08)

##############################################################################
#
# Fixation probability panels
#
##############################################################################

frequency_xticks = [0,0.2,0.4,0.6,0.8,1]

pfix_axis = plt.Subplot(fig, inner_grid_1[0])
fig.add_subplot(pfix_axis)

pfix_axis.set_ylabel('Survival probability')
#pfix_axis.set_xlabel('Observed frequency')

pfix_axis.spines['top'].set_visible(False)
pfix_axis.spines['right'].set_visible(False)
pfix_axis.get_xaxis().tick_bottom()
pfix_axis.get_yaxis().tick_left()
pfix_axis.get_yaxis().set_tick_params(direction='out')
pfix_axis.get_xaxis().set_tick_params(direction='out')

pfix_axis.set_xticks(frequency_xticks)
pfix_axis.set_xticklabels([])

pfix_axis.set_xlim([0,1.01])
pfix_axis.set_ylim([0,1.01])
line, = pfix_axis.plot([0,1],[1,1],'k:',linewidth=0.25, label='Hitchhiking')
line.set_dashes((0.5,0.5))
line, = pfix_axis.plot([0,1],[0,1],'k--',linewidth=0.25, label='Quasi-neutral')
line.set_dashes((2,1))

pooled_pfix_axis = plt.Subplot(fig, inner_grid_1[1])
fig.add_subplot(pooled_pfix_axis)

pooled_pfix_axis.set_ylabel('Survival probability')
pooled_pfix_axis.set_xlabel('Allele frequency (within clade)')

pooled_pfix_axis.spines['top'].set_visible(False)
pooled_pfix_axis.spines['right'].set_visible(False)
pooled_pfix_axis.get_xaxis().tick_bottom()
pooled_pfix_axis.get_yaxis().tick_left()
pooled_pfix_axis.get_yaxis().set_tick_params(direction='out')
pooled_pfix_axis.get_xaxis().set_tick_params(direction='out')

pooled_pfix_axis.set_xticks(frequency_xticks)
pooled_pfix_axis.set_xlim([0,1.05])
pooled_pfix_axis.set_ylim([0,1.05])
line, = pooled_pfix_axis.plot([0,1],[1,1],'k:',linewidth=0.25)
line.set_dashes((0.5, 0.5))
line, = pooled_pfix_axis.plot([0,1],[0,1],'k--',linewidth=0.25)
line.set_dashes((2,1))

##############################################################################
#
# Fixation trajectory panel
#
##############################################################################

fixation_trajectory_axis = plt.Subplot(fig, inner_grid_2[0])
fig.add_subplot(fixation_trajectory_axis)

fixation_trajectory_axis.set_ylabel('Fixed mutations (major clade)')
#fixation_trajectory_axis.set_xlabel('Generation')

fixation_trajectory_axis.spines['top'].set_visible(False)
fixation_trajectory_axis.spines['right'].set_visible(False)
fixation_trajectory_axis.get_xaxis().tick_bottom()
fixation_trajectory_axis.get_yaxis().tick_left()
fixation_trajectory_axis.get_yaxis().set_tick_params(direction='out')
fixation_trajectory_axis.get_xaxis().set_tick_params(direction='out')

fixation_trajectory_axis.set_xlim([-1000,61000])
fixation_trajectory_axis.set_ylim([0,110])


xticks = [10000*i for i in xrange(0,7)]
xticklabels = ["" for i in xrange(0,7)]
fixation_trajectory_axis.set_xticks(xticks)
fixation_trajectory_axis.set_xticklabels(xticklabels)

##############################################################################
#
# Transit time panel
#
##############################################################################


fixation_time_axis = plt.Subplot(fig, inner_grid_2[1])
fig.add_subplot(fixation_time_axis)

fixation_time_axis.set_ylabel('Transit time')
fixation_time_axis.set_xlabel('Generation')

#fixation_time_axis.spines['top'].set_visible(False)
fixation_time_axis.spines['right'].set_visible(False)
fixation_time_axis.get_xaxis().tick_bottom()
fixation_time_axis.get_yaxis().tick_left()
fixation_time_axis.get_yaxis().set_tick_params(direction='out')
fixation_time_axis.get_xaxis().set_tick_params(direction='out')

fixation_time_axis.set_xlim([-1000,61000])
fixation_time_axis.set_ylim([0,10000])

yticks = [1000*2*i for i in xrange(0,6)]
yticklabels = ['%dk' % (i*2) for i in xrange(0,6)]
yticklabels[0] = '0'
fixation_time_axis.set_yticks(yticks)
fixation_time_axis.set_yticklabels(yticklabels)

xticks = [10000*i for i in xrange(0,7)]
xticklabels = ['%dk' % (10*i) for i in xrange(0,7)]
xticklabels[0] = '0'
fixation_time_axis.set_xticks(xticks)
fixation_time_axis.set_xticklabels(xticklabels)

######
#
# Now do the calculations
#
######

#populations = ['m5','p2','p4','p1','p5','p6','p3','m4','m2','m1','m3','m6']
nonmutator_populations = ['m5','p2','p4','p1','m6','p5']
mutator_populations = ['m1','m2','m3','m4','p3','p6']

population_labels = {population: parse_file.get_pretty_name(population) for population in nonmutator_populations}
population_colors = { nonmutator_populations[i] : parse_file.nonmutator_line_colors[i] for i in xrange(0,len(nonmutator_populations)) } 

for i in xrange(0,len(mutator_populations)):
    population_labels[mutator_populations[i]] = parse_file.get_pretty_name(mutator_populations[i])
    population_colors[mutator_populations[i]] = parse_file.mutator_line_colors[i]

populations = mutator_populations+nonmutator_populations

frequencies = numpy.linspace(0.1,0.9,201)
df = 0.05
fstar = 0.5

tstar = 20025

nonmutator_num_in_bin = numpy.zeros_like(frequencies)
nonmutator_avg_f = numpy.zeros_like(frequencies)

late_nonmutator_num_in_bin = numpy.zeros_like(frequencies)
late_nonmutator_avg_f = numpy.zeros_like(frequencies)

mutator_num_in_bin = numpy.zeros_like(frequencies)
mutator_avg_f = numpy.zeros_like(frequencies)

late_mutator_num_in_bin = numpy.zeros_like(frequencies)
late_mutator_avg_f = numpy.zeros_like(frequencies)

origin_fixation_times = {population: ([],[],[]) for population in populations}


for population in populations:

    sys.stderr.write("Processing fixation probabilities for %s...\t" % population)

    is_mutator = (population in parse_file.mutator_lines)

    num_in_bin = numpy.zeros_like(nonmutator_num_in_bin)
    avg_f = numpy.zeros_like(nonmutator_avg_f)

    # load mutations
    mutations, depth_tuple = parse_file.parse_annotated_timecourse(population)
    population_avg_depth_times, population_avg_depths, clone_avg_depth_times, clone_avg_depths = depth_tuple
    dummy_times,fmajors,fminors,haplotype_trajectories = parse_file.parse_haplotype_timecourse(population)
    state_times, state_trajectories = parse_file.parse_well_mixed_state_timecourse(population)
    
    fremainders = 1-fmajors-fminors
    
    ffixeds = numpy.ones_like(fmajors)*1.0
    
    num_runs = []
    
    # we want to be conservative, so use the largest possible
    # parent clade to estimate within-clade freq
    fextincts = numpy.fmax(fremainders,numpy.fmax(fmajors,fminors))
    
    for mutation_idx in xrange(0,len(mutations)):

        location, gene_name, allele, var_type, test_statistic, pvalue, cutoff_idx, depth_fold_change, depth_change_pvalue, times, alts, depths, clone_times, clone_alts, clone_depths = mutations[mutation_idx] 

        if is_mutator and var_type=='sv' or var_type=='indel':
            continue
    
        Ls = haplotype_trajectories[mutation_idx]
        state_Ls = state_trajectories[mutation_idx]
        
        good_idxs, filtered_alts, filtered_depths = timecourse_utils.mask_timepoints(times, alts, depths, var_type, cutoff_idx, depth_fold_change, depth_change_pvalue)
    
        freqs = timecourse_utils.estimate_frequencies(filtered_alts, filtered_depths)
        
        masked_times = times[good_idxs]
        masked_freqs = freqs[good_idxs]
        masked_depths = depths[good_idxs]
        masked_state_Ls = state_Ls[good_idxs]
        masked_Ls = Ls[good_idxs]
        
        masked_fmajors = fmajors[good_idxs]
        masked_fminors = fminors[good_idxs]
        masked_ffixeds = ffixeds[good_idxs]
        masked_fextincts = fextincts[good_idxs]
        
        majority_fixed_states = set([parse_file.clade_hmm_states['FB'], parse_file.clade_hmm_states['FM'], parse_file.clade_hmm_states['PB*']])

        # Estimate appearance and fixation times      
        if masked_Ls[-1] in majority_fixed_states:
            t0,tf,dt = timecourse_utils.calculate_appearance_fixation_time_from_clade_hmm(masked_times, masked_freqs, masked_Ls)
            origin_fixation_times[population][0].append(t0)
            origin_fixation_times[population][1].append(tf)
            origin_fixation_times[population][2].append(dt)
                
        # Now split the trajectory into independent polymorphic runs
        # each of which contains a single final state (unless end of experiment)
        independent_runs = timecourse_utils.split_clade_hmm(masked_times,masked_freqs, masked_Ls)
        
        num_runs.append(len(independent_runs))
            
        for run_idxs in independent_runs:
            
            if len(run_idxs)<2:
                # need at least one polymorphic state and one final state
                continue 
            
            # initial time
            t = masked_times[run_idxs[0]]
            
            # get final state
            final_state = masked_Ls[run_idxs[-1]]
            
            # get frequency of parent clade during run
            if final_state == parse_file.clade_hmm_states['FB'] or final_state==parse_file.clade_hmm_states['PB'] or final_state==parse_file.clade_hmm_states['PB*']:
                parent_freqs = masked_ffixeds
            elif final_state == parse_file.clade_hmm_states['FM'] or final_state == parse_file.clade_hmm_states['PM']:
                parent_freqs = masked_fmajors
            elif final_state == parse_file.clade_hmm_states['Fm'] or final_state == parse_file.clade_hmm_states['Pm']:
                parent_freqs = masked_fmajors      
            else:
                parent_freqs = masked_fextincts     
      
            renormalized_freqs = numpy.clip(masked_freqs[run_idxs]/parent_freqs[run_idxs],0,1)
            
            # get fixation weight
            if final_state in parse_file.clade_fixed_states:
                fixation_weight = 1.0
            elif final_state in parse_file.clade_extinct_states:
                fixation_weight = 0.0
            else:
                fixation_weight = renormalized_freqs[-1] > fstar
            
            #individual_bin_weights = (renormalized_freqs[:,None]>=(frequencies[None,:]-df/2))*(renormalized_freqs[:,None]<=(frequencies[None,:]+df/2))
            
            individual_bin_weights = numpy.exp(-numpy.power((renormalized_freqs[:,None]-frequencies[None,:])/df,2))
            individual_bin_weights *= (individual_bin_weights>0.01)
            
            bin_weights = individual_bin_weights.sum(axis=0)
            fixation_weights = (individual_bin_weights*fixation_weight).sum(axis=0)
              
            num_in_bin += bin_weights
            avg_f += fixation_weights
            
            if population in parse_file.complete_nonmutator_lines: 
                
                nonmutator_num_in_bin += bin_weights
                nonmutator_avg_f += fixation_weights
                
                if t > tstar:
                    late_nonmutator_num_in_bin += bin_weights
                    late_nonmutator_avg_f += fixation_weights
            else:
                
                mutator_num_in_bin += bin_weights
                mutator_avg_f += fixation_weights
                
                if t > tstar:
                    late_mutator_num_in_bin += bin_weights
                    late_mutator_avg_f += fixation_weights
             
    
    avg_f = avg_f/(num_in_bin+(num_in_bin<1))
        
    pfix_axis.plot(frequencies[(num_in_bin>=1)], avg_f[(num_in_bin>=1)], '-', color=parse_file.get_line_color(population),linewidth=0.5,markersize=1)

    sys.stderr.write("Done!\n")
    
nonmutator_avg_f = nonmutator_avg_f/nonmutator_num_in_bin
late_nonmutator_avg_f = late_nonmutator_avg_f/late_nonmutator_num_in_bin
mutator_avg_f = mutator_avg_f/mutator_num_in_bin
late_mutator_avg_f = late_mutator_avg_f/late_mutator_num_in_bin


pooled_pfix_axis.plot(frequencies, late_mutator_avg_f,'-',color=mutator_color,alpha=0.5,linewidth=0.5,markersize=1)

pooled_pfix_axis.plot(frequencies, mutator_avg_f,'-',color=mutator_color,linewidth=0.5,markersize=1)

pooled_pfix_axis.plot(frequencies, late_nonmutator_avg_f,'-',color=nonmutator_color,linewidth=0.5,alpha=0.5,markersize=1)

pooled_pfix_axis.plot(frequencies, nonmutator_avg_f,'-',color=nonmutator_color,linewidth=0.5,markersize=1)

# set up legends (offscreen)

pooled_pfix_axis.plot([-2,-1],[-2,-1], '-', color=nonmutator_color, label='Wildtype (all)',markersize=1,markeredgewidth=0,linewidth=0.5)

pooled_pfix_axis.plot([-2,-1],[-2,-1] ,'-', color=nonmutator_color, label='Wildtype (>20k)',markersize=1,markeredgewidth=0,alpha=0.5,linewidth=0.5)

pooled_pfix_axis.plot([-2,-1],[-2,-1],'-', color=mutator_color, label='Mutator (all)',markersize=1,markeredgewidth=0,linewidth=0.5)

pooled_pfix_axis.plot([-2,-1],[-2,-1],'-', color=mutator_color, markersize=1,markeredgewidth=0,alpha=0.5,linewidth=0.5, label='Mutator (>20k)')

pfix_axis.legend(loc='lower right',frameon=False,fontsize=4,numpoints=1)   
pooled_pfix_axis.legend(loc='lower right',frameon=False,fontsize=4,numpoints=1)   


############
#
# Now do within-clade transit time calculation
#
############


masked_populations = set([])

total_origin_times = []
total_fixation_intervals = []

for population in parse_file.complete_nonmutator_lines:
    
    sys.stderr.write('Processing %s...\n' % population)
    
    color = parse_file.get_line_color(population)
    
    if population in origin_fixation_times.keys() and (population not in masked_populations):
    
        origin_times = numpy.array(origin_fixation_times[population][0])
        fixation_times = numpy.array(origin_fixation_times[population][1])
        fixation_intervals = numpy.array(origin_fixation_times[population][2])
        
        total_origin_times.extend(origin_times)
        total_fixation_intervals.extend(fixation_intervals)
            
        fixation_time_axis.plot(origin_times, fixation_intervals, '.', label=population,color=color,markersize=2)

        
total_origin_times, total_fixation_intervals = (numpy.array(x) for x in zip(*sorted(zip(total_origin_times, total_fixation_intervals), key=lambda pair: (pair[0]))))    

# sliding window medians
mid_windows = numpy.linspace(0,60000,50)
lower_windows = mid_windows-5000
upper_windows = mid_windows+5000

avg_dts = []
upper_dts = []
lower_dts = []
avg_times = []

for tlower,tmid,tupper in zip(lower_windows,mid_windows,upper_windows):
    
    fixation_intervals = total_fixation_intervals[(total_origin_times>=tlower)*(total_origin_times<=tupper)]
    if len(fixation_intervals)>0:
    
    
        avg_times.append(tmid)
        avg_dts.append(numpy.median(fixation_intervals))
    
        upper_dts.append( numpy.percentile(fixation_intervals,75,interpolation='nearest')) 
        lower_dts.append( numpy.percentile(fixation_intervals,25,interpolation='nearest'))
    

fixation_time_axis.plot(avg_times, avg_dts, 'w-', linewidth=1,path_effects=[pe.Stroke(linewidth=2.5, foreground='k'), pe.Normal()])


fixation_time_axis.fill_between(avg_times, lower_dts, upper_dts,color='0.8')
fixation_time_axis.plot(avg_times, lower_dts,'k-',linewidth=0.25)
fixation_time_axis.plot(avg_times, upper_dts,'k-',linewidth=0.25)  

######
#
# Now do fixation trajectories
#
########

theory_times = numpy.arange(0,121)*500
total_Mfixeds = numpy.zeros_like(theory_times)*1.0

num_pops = 0.0
for population in parse_file.complete_nonmutator_lines:
    
    color=parse_file.get_line_color(population)
    
    if population in origin_fixation_times.keys():
        num_pops+=1
        # construct timecourse
        fixation_times = numpy.array(origin_fixation_times[population][1])
    
        Mfixeds = numpy.zeros_like(total_Mfixeds)
        for tf in fixation_times:
            Mfixeds += (theory_times>=tf)
            
            #idx = long((t+250)/500)
            #if idx>len(theory_times):
            #    continue
            #Mfixeds[idx:] += 1

        total_Mfixeds += Mfixeds
        
        fixation_trajectory_axis.plot(theory_times, Mfixeds,color=color)

avg_Mfixeds = total_Mfixeds/num_pops

fixation_trajectory_axis.plot(theory_times, avg_Mfixeds, 'w-' ,linewidth=1,path_effects=[pe.Stroke(linewidth=2.5, foreground='k'), pe.Normal()])

############
#
# Add figure labels
#
############

fixation_trajectory_axis.text(1500, 95, figure_utils.get_panel_label('a'),fontsize=6,fontweight='bold')
fixation_time_axis.text(1500, 9000, figure_utils.get_panel_label('b'), fontsize=6, fontweight='bold')
pfix_axis.text(0.05, 0.9, figure_utils.get_panel_label('c'), fontsize=6, fontweight='bold')
pooled_pfix_axis.text(0.05, 0.9, figure_utils.get_panel_label('d'), fontsize=6, fontweight='bold')

############
#
# Save figure
#
############
sys.stderr.write("Saving figure...\t")
fig.savefig(parse_file.figure_directory+'fig4.pdf', bbox_inches='tight')
sys.stderr.write("Done!\n")