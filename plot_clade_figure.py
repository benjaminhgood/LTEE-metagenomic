###############################
#
# Rest of script begins here
#
################################

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
import bz2
import parse_file
import matplotlib
import matplotlib.pyplot as plt
import timecourse_utils
import figure_utils

focal_population = 'm6'
#focal_population = 'p5'
#focal_population = 'm4'
#focal_population = 'm5'
#focal_population = 'm1'
all_populations = parse_file.complete_nonmutator_lines + parse_file.mutator_lines
all_colors = parse_file.nonmutator_line_colors+parse_file.mutator_line_colors

remaining_populations = []
remaining_colors = []
for population,color in zip(all_populations, all_colors):
    if population!=focal_population:
        remaining_populations.append(population)
        remaining_colors.append(color)

remaining_populations = all_populations
remaining_colors = all_colors

##############################
#
# First set up the figure
#
##############################

mpl.rcParams['font.size'] = 4.0
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'

generation_ticks = [10000*i for i in xrange(0,7)]
generation_ticklabels = ['%dk' % (10*i) for i in xrange(0,7)]
generation_ticklabels[0] = '0'

generation_ticks = [5000*i for i in xrange(0,13)]
generation_ticklabels = ['%dk' % (5*i) for i in xrange(0,13)]

frequency_ticks = [0,0.2,0.4,0.6,0.8,1.0]

pylab.figure(1,figsize=(7.2,3.5))
fig = pylab.gcf()

# make three panels panels
outer_grid  = gridspec.GridSpec(1, 2, width_ratios=[4,1], wspace=0.15)

inner_grid_1 = gridspec.GridSpecFromSubplotSpec(4, 1, height_ratios=[1,1,1,1],
                subplot_spec=outer_grid[0], hspace=0.25) #, hspace=0.08)

inner_grid_2 = gridspec.GridSpecFromSubplotSpec(len(remaining_populations), 1, height_ratios=([1]*len(remaining_populations)),
                subplot_spec=outer_grid[1], hspace=0.15) #, hspace=0.08)

fixed_axis = plt.Subplot(fig, inner_grid_1[0])
fig.add_subplot(fixed_axis)

fixed_axis.set_title('%s basal clade' % parse_file.get_pretty_name(focal_population),loc='right',fontsize=4,y=0.91)
    

fixed_axis.spines['top'].set_visible(False)
fixed_axis.spines['right'].set_visible(False)
fixed_axis.get_xaxis().tick_bottom()
fixed_axis.get_yaxis().tick_left()
fixed_axis.get_yaxis().set_tick_params(direction='out',length=3,pad=1)
fixed_axis.get_xaxis().set_tick_params(direction='out',length=3,pad=1)

fixed_axis.set_xticks(generation_ticks)
fixed_axis.set_xticklabels([])
fixed_axis.set_yticks(frequency_ticks)

fixed_axis.set_ylabel('Allele frequency, $f(t)$')


fixed_axis.set_ylim([0,1.01])
fixed_axis.set_xlim([0,60100])

majority_axis = plt.Subplot(fig, inner_grid_1[1])
fig.add_subplot(majority_axis)

majority_axis.set_title('%s major clade' % parse_file.get_pretty_name(focal_population),loc='right',fontsize=4,y=0.91)


majority_axis.spines['top'].set_visible(False)
majority_axis.spines['right'].set_visible(False)
majority_axis.get_xaxis().tick_bottom()
majority_axis.get_yaxis().tick_left()
majority_axis.get_yaxis().set_tick_params(direction='out',length=3,pad=1)
majority_axis.get_xaxis().set_tick_params(direction='out',length=3,pad=1)

majority_axis.set_xticks(generation_ticks)
majority_axis.set_xticklabels([])
majority_axis.set_yticks(frequency_ticks)
majority_axis.set_ylabel('Allele frequency, $f(t)$')

majority_axis.set_ylim([0,1.01])
majority_axis.set_xlim([0,60100])


minority_axis = plt.Subplot(fig, inner_grid_1[2])
fig.add_subplot(minority_axis)

minority_axis.set_title('%s minor clade' % parse_file.get_pretty_name(focal_population),loc='right',fontsize=4,y=0.91)


minority_axis.spines['top'].set_visible(False)
minority_axis.spines['right'].set_visible(False)
minority_axis.get_xaxis().tick_bottom()
minority_axis.get_yaxis().tick_left()
minority_axis.get_yaxis().set_tick_params(direction='out',length=3,pad=1)
minority_axis.get_xaxis().set_tick_params(direction='out',length=3,pad=1)

minority_axis.set_xticks(generation_ticks)
minority_axis.set_xticklabels([])
minority_axis.set_yticks(frequency_ticks)

minority_axis.set_ylabel('Allele frequency, $f(t)$')


minority_axis.set_ylim([0,1.01])
minority_axis.set_xlim([0,60100])


extinct_axis = plt.Subplot(fig, inner_grid_1[3])
fig.add_subplot(extinct_axis)

extinct_axis.set_title('%s extinct' % parse_file.get_pretty_name(focal_population),loc='right',fontsize=4,y=0.91)


extinct_axis.set_xlabel('Generation, $t$')

extinct_axis.spines['top'].set_visible(False)
extinct_axis.spines['right'].set_visible(False)
extinct_axis.get_xaxis().tick_bottom()
extinct_axis.get_yaxis().tick_left()
extinct_axis.get_yaxis().set_tick_params(direction='out',length=3,pad=1)
extinct_axis.get_xaxis().set_tick_params(direction='out',length=3,pad=1)

extinct_axis.set_xticks(generation_ticks)
extinct_axis.set_xticklabels(generation_ticklabels)
extinct_axis.set_yticks(frequency_ticks)

extinct_axis.set_ylabel('Allele frequency, $f(t)$')


extinct_axis.set_ylim([0,1.01])
extinct_axis.set_xlim([0,60100])

remaining_axes = []
for i in xrange(0, len(remaining_populations)):
    remaining_axis = plt.Subplot(fig, inner_grid_2[i])
    fig.add_subplot(remaining_axis)
    remaining_axes.append(remaining_axis)
    remaining_axis.set_ylabel(parse_file.get_pretty_name(remaining_populations[i]),fontsize=4,labelpad=-2)
    
    remaining_axis.spines['top'].set_visible(False)
    remaining_axis.spines['right'].set_visible(False)
    remaining_axis.get_xaxis().tick_bottom()
    remaining_axis.get_yaxis().tick_left()
    remaining_axis.get_yaxis().set_tick_params(direction='out',length=2)
    remaining_axis.get_xaxis().set_tick_params(direction='out',length=2)
    
    remaining_axis.set_ylim([0,1])
    remaining_axis.set_xlim([0,61000])
    
    if i==(len(remaining_populations)-1):
        remaining_axis.set_xticks([i*10000 for i in xrange(0,7)])
        remaining_axis.set_xticklabels([])
    else:
        remaining_axis.set_xticks([])
    
    remaining_axis.set_yticks([0,0.5,1])
    remaining_axis.set_yticklabels([])

########################################
#
# Now do the plotting (focal first, then rest)
#
########################################


        
theory_times = numpy.arange(0,121)*500

gene_names, start_positions, end_positions, promoter_start_positions, promoter_end_positions, gene_sequences, strands = parse_file.parse_gene_list()

gene_name_position_map = {gene_names[i]: (start_positions[i], end_positions[i]) for i in xrange(0,len(gene_names))}

sys.stderr.write('Processing focal population: %s...\t' % focal_population)

# load mutations
mutations, depth_tuple = parse_file.parse_annotated_timecourse(focal_population)
    
population_avg_depth_times, population_avg_depths, clone_avg_depth_times, clone_avg_depths = depth_tuple
    
dummy_times,fmajors,fminors,haplotype_trajectories = parse_file.parse_haplotype_timecourse(focal_population)
   
is_mutator = (focal_population in parse_file.mutator_lines)
    
for mutation_idx in xrange(0,len(mutations)):
 
    location, gene_name, allele, var_type, test_statistic, pvalue, cutoff_idx, depth_fold_change, depth_change_pvalue, times, alts, depths, clone_times, clone_alts, clone_depths = mutations[mutation_idx] 
    
    if is_mutator and var_type=='sv' or var_type=='indel':
        continue
    
    Ls = haplotype_trajectories[mutation_idx]
    
    good_idxs, filtered_alts, filtered_depths = timecourse_utils.mask_timepoints(times, alts, depths, var_type, cutoff_idx, depth_fold_change, depth_change_pvalue)
        
    freqs = timecourse_utils.estimate_frequencies(filtered_alts, filtered_depths)
        
    masked_times = times[good_idxs]
    masked_freqs = freqs[good_idxs]
        
    interpolation_function = timecourse_utils.create_interpolation_function(masked_times, masked_freqs, tmax=100000)
    haplotype_interpolation_function = timecourse_utils.create_interpolation_function(times, Ls, tmax=100000,kind='nearest')
    
    clone_freqs = clone_alts*1.0/(clone_depths+(clone_depths==0))
    clone_depth_fold_changes = timecourse_utils.estimate_depth_fold_changes(clone_avg_depths, clone_depths)

    # classify clade
    if (Ls[-1]==parse_file.ANCESTRAL_FIXED) or (Ls[-1]==parse_file.ANCESTRAL_POLYMORPHIC):
        clade = "ancestral"
    elif (Ls[-1]==parse_file.MAJOR_FIXED) or (Ls[-1]==parse_file.MAJOR_POLYMORPHIC):
        clade = "majority"
    elif (Ls[-1]==parse_file.MINOR_FIXED) or (Ls[-1]==parse_file.MINOR_POLYMORPHIC):
        clade = "minority"
    else:
        clade = "extinct"
    
        
    
      
    if clade=='ancestral':
        fixed_axis.plot(masked_times, masked_freqs, '-o', alpha=0.5, markersize=1, markeredgecolor='none', zorder=4, linewidth=0.5)
    else:  
        fixed_axis.plot(masked_times, masked_freqs, '-', color='0.7', alpha=0.5, linewidth=0.25,zorder=1)
    
    if clade=='majority':
        majority_axis.plot(masked_times, masked_freqs, '-o', alpha=0.5, markersize=1, markeredgecolor='none', zorder=4, linewidth=0.5)
    else:  
        majority_axis.plot(masked_times, masked_freqs, '-', color='0.7', alpha=0.5, linewidth=0.25,zorder=1)
    
    if clade=='minority':
        minority_axis.plot(masked_times, masked_freqs, '-o', alpha=0.5, markersize=1, markeredgecolor='none', zorder=4, linewidth=0.5)
    else:  
        minority_axis.plot(masked_times, masked_freqs, '-', color='0.7', alpha=0.5, linewidth=0.25,zorder=1)
    
    if clade=='extinct':
        extinct_axis.plot(masked_times, masked_freqs, '-o', alpha=0.5, markersize=1, markeredgecolor='none', zorder=4, linewidth=0.5)
    else:  
        extinct_axis.plot(masked_times, masked_freqs, '-', color='0.7', alpha=0.5, linewidth=0.25,zorder=1)
    
majority_axis.plot(dummy_times, fmajors, 'k-', linewidth=2,zorder=10)
minority_axis.plot(dummy_times, fminors, 'k-', linewidth=2,zorder=10)

sys.stderr.write('Done!\n')

for idx in xrange(0,len(remaining_populations)):
    population = remaining_populations[idx]
    freq_axis = remaining_axes[idx]
    color = remaining_colors[idx]
    
    sys.stderr.write('Processing sidebar population %s...\t' % population)

    # load mutations
    mutations, depth_tuple = parse_file.parse_annotated_timecourse(population)
    
    population_avg_depth_times, population_avg_depths, clone_avg_depth_times, clone_avg_depths = depth_tuple
    
    dummy_times,fmajors,fminors,haplotype_trajectories = parse_file.parse_haplotype_timecourse(population)
    
    for mutation_idx in xrange(0,len(mutations)):
 
        location, gene_name, allele, var_type, test_statistic, pvalue, cutoff_idx, depth_fold_change, depth_change_pvalue, times, alts, depths, clone_times, clone_alts, clone_depths = mutations[mutation_idx] 
    
        Ls = haplotype_trajectories[mutation_idx]
    
        good_idxs, filtered_alts, filtered_depths = timecourse_utils.mask_timepoints(times, alts, depths, var_type, cutoff_idx, depth_fold_change, depth_change_pvalue)
        
        freqs = timecourse_utils.estimate_frequencies(filtered_alts, filtered_depths)
        
        masked_times = times[good_idxs]
        masked_freqs = freqs[good_idxs]
        
        interpolation_function = timecourse_utils.create_interpolation_function(masked_times, masked_freqs, tmax=100000)
        haplotype_interpolation_function = timecourse_utils.create_interpolation_function(times, Ls, tmax=100000,kind='nearest')
    
        clone_freqs = clone_alts*1.0/(clone_depths+(clone_depths==0))
        clone_depth_fold_changes = timecourse_utils.estimate_depth_fold_changes(clone_avg_depths, clone_depths)
    
        freq_axis.plot(masked_times, masked_freqs, '-', color='0.7', alpha=0.5,linewidth=0.25,zorder=1)
        #freq_axis.plot(masked_times, masked_freqs, '-', alpha=0.5,linewidth=0.25,zorder=1)

    freq_axis.plot(dummy_times, fminors, '-', color='#f768a1', linewidth=1,zorder=10)     
    freq_axis.plot(dummy_times, fmajors, '-', color='#7a0177', linewidth=1,zorder=10)

    sys.stderr.write("Done!\n")

#adding figure labels
fixed_axis.text(1000, 1.05, figure_utils.get_panel_label('a'),fontsize=6,fontweight='bold')
remaining_axes[0].text(0, 1.15, figure_utils.get_panel_label('b'),fontsize=6,fontweight='bold')

sys.stderr.write("Saving final PNG image...\t")
fig.savefig(parse_file.figure_directory+'fig3.png',bbox_inches='tight', dpi=300, transparent=True)
pylab.close(fig)
sys.stderr.write("Done!\n")