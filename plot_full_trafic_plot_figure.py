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

debug=False
remove_data_for_copyediting = False

focal_populations = ['p2','m6','m1']
all_populations = parse_file.complete_nonmutator_lines + parse_file.mutator_lines
all_colors = parse_file.nonmutator_line_colors+parse_file.mutator_line_colors

remaining_populations = ['m5','p1','p4','p5','m2','m3','m4','p3','p6']

##############################
#
# First set up the figure
#
##############################

mpl.rcParams['font.size'] = 5.0
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'
mpl.rcParams['axes.labelpad'] = 2
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['font.sans-serif'] = 'Arial'

mpl.rcParams['font.serif'] = 'Times New Roman'
mpl.rcParams['mathtext.rm'] = 'serif'
mpl.rcParams['mathtext.it'] = 'serif:italic'
mpl.rcParams['mathtext.bf'] = 'serif:bold'
mpl.rcParams['mathtext.fontset'] = 'custom'



generation_ticks = [10000*i for i in xrange(0,7)]
generation_ticklabels = ['%dk' % (10*i) for i in xrange(0,7)]
generation_ticklabels[0] = '0'

generation_ticks = [5000*i for i in xrange(0,13)]
generation_ticklabels = ['%dk' % (5*i) for i in xrange(0,13)]

frequency_ticks = [0,0.2,0.4,0.6,0.8,1.0]

#pylab.figure(1,figsize=(7.09,2.7))
pylab.figure(1,figsize=(8.5,2.9))
fig = pylab.gcf()

# make three panels panels
outer_grid  = gridspec.GridSpec(1, 2, width_ratios=[4,1], wspace=0.15)

inner_grid_1 = gridspec.GridSpecFromSubplotSpec(len(focal_populations), 1, height_ratios=([1]*len(focal_populations)),
                subplot_spec=outer_grid[0], hspace=0.25) #, hspace=0.08)

inner_grid_2 = gridspec.GridSpecFromSubplotSpec(len(remaining_populations), 1, height_ratios=([1]*len(remaining_populations)),
                subplot_spec=outer_grid[1], hspace=0.2) #, hspace=0.08)


focal_axes = []
for i in xrange(0,len(focal_populations)):
    axis = plt.Subplot(fig, inner_grid_1[i])
    focal_axes.append(axis)
    fig.add_subplot(axis)

    axis.set_title('%s' % parse_file.get_pretty_name(focal_populations[i]),loc='right',fontsize=5,y=0.93)
    

    axis.spines['top'].set_visible(False)
    axis.spines['right'].set_visible(False)
    axis.get_xaxis().tick_bottom()
    axis.get_yaxis().tick_left()
    axis.get_yaxis().set_tick_params(direction='out',length=3,pad=1)
    axis.get_xaxis().set_tick_params(direction='out',length=3,pad=1)   

    axis.set_xticks(generation_ticks)
    axis.set_yticks(frequency_ticks)
    
    if i==(len(focal_populations)-1):
        axis.set_xlabel('Generation, $t$')
        axis.set_xticklabels(generation_ticklabels)
    else:
        axis.set_xticklabels([])
        
    axis.set_ylim([0,1.01])
    axis.set_xlim([0,60100])

    axis.set_ylabel('Allele frequency, $f(t)$')


#focal_axes[1].set_ylabel('Allele frequency, $f(t)$')

remaining_axes = []
for i in xrange(0, len(remaining_populations)):
    remaining_axis = plt.Subplot(fig, inner_grid_2[i])
    fig.add_subplot(remaining_axis)
    remaining_axes.append(remaining_axis)
    remaining_axis.set_ylabel(parse_file.get_pretty_name(remaining_populations[i]),fontsize=5,labelpad=-2)
    
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

if debug==True:
    starting_idx = -5
else:
    starting_idx = 0
   
theory_times = numpy.arange(0,121)*500

gene_names, start_positions, end_positions, promoter_start_positions, promoter_end_positions, gene_sequences, strands = parse_file.parse_gene_list()

gene_name_position_map = {gene_names[i]: (start_positions[i], end_positions[i]) for i in xrange(0,len(gene_names))}
  
for i in xrange(0,len(focal_populations)):

    sys.stderr.write('Processing focal population: %s...\t' % focal_populations[i])

    # load mutations
    mutations, depth_tuple = parse_file.parse_annotated_timecourse(focal_populations[i])
    
    population_avg_depth_times, population_avg_depths, clone_avg_depth_times, clone_avg_depths = depth_tuple
    
    dummy_times,fmajors,fminors,haplotype_trajectories = parse_file.parse_haplotype_timecourse(focal_populations[i])
    
    for mutation_idx in range(0,len(mutations))[starting_idx:]:
 
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
        
        if not remove_data_for_copyediting:
            focal_axes[i].plot(masked_times, masked_freqs, '-o', alpha=0.5, markersize=1, markeredgecolor='none', zorder=4, linewidth=0.5) #,rasterized=True)
    
    sys.stderr.write("Done!\n")
    
for idx in xrange(0,len(remaining_populations)):
    population = remaining_populations[idx]
    freq_axis = remaining_axes[idx]
    #color = remaining_colors[idx]
    
    if False: #population in focal_populations:
        freq_axis.spines['bottom'].set_visible(False)
        freq_axis.spines['left'].set_visible(False)
        freq_axis.set_xticks([])
        freq_axis.set_yticks([])
        freq_axis.set_ylabel("")
        continue
    
    sys.stderr.write('Processing sidebar population %s...\t' % population)

    # load mutations
    mutations, depth_tuple = parse_file.parse_annotated_timecourse(population)
    
    population_avg_depth_times, population_avg_depths, clone_avg_depth_times, clone_avg_depths = depth_tuple
    
    dummy_times,fmajors,fminors,haplotype_trajectories = parse_file.parse_haplotype_timecourse(population)
    
    for mutation_idx in range(0,len(mutations))[starting_idx:]:
 
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
    
        
        if not remove_data_for_copyediting:
            freq_axis.plot(masked_times, masked_freqs, '-', alpha=0.5,linewidth=0.25,zorder=1) #,rasterized=True)

    sys.stderr.write("Done!\n")
    


if remove_data_for_copyediting:
    figure_name = 'fig1_nodata'
else:
    figure_name = 'fig1'
    
sys.stderr.write("Saving figure...\t")
fig.savefig(parse_file.figure_directory+figure_name+'.png',bbox_inches='tight', dpi=600)
fig.savefig(parse_file.figure_directory+figure_name+'.pdf',bbox_inches='tight') #, dpi=300)
pylab.close(fig)
sys.stderr.write("Done!\n")
sys.exit(1)
