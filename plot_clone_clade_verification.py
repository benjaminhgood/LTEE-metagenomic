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

# dictionary between population and dictionary between clone and list of classification of mutations that it "has" at that timepoint. do not add anything if it is deleted. 

# clone_clades[population][clone] 
# then figure is bar charts (sorted by timepoints) 
# separate line for each population

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


metapopulations = [parse_file.complete_nonmutator_lines, parse_file.mutator_lines]


##############################
#
# First set up the figure
#
##############################

mpl.rcParams['font.size'] = 5.0
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'

pylab.figure(1,figsize=(7,6.5))
fig = pylab.gcf()

# make three panels panels
outer_grid  = gridspec.GridSpec(1, 2, wspace=0.25)

inner_grid_1 = gridspec.GridSpecFromSubplotSpec(6, 1, height_ratios=([1]*6),
                subplot_spec=outer_grid[0], hspace=0.2) #, hspace=0.08)

inner_grid_2 = gridspec.GridSpecFromSubplotSpec(6, 1, height_ratios=([1]*6),
                subplot_spec=outer_grid[1], hspace=0.2) #, hspace=0.08)

inner_grids = [inner_grid_1, inner_grid_2]

axes = []
for metapopulation_idx in xrange(0,2):
    axes.append([])
    for population_idx in xrange(0,6):
        
        population = metapopulations[metapopulation_idx][population_idx]
        
        axis = plt.Subplot(fig, inner_grids[metapopulation_idx][population_idx])
        fig.add_subplot(axis)
        axes[metapopulation_idx].append(axis)
        axis.set_title( parse_file.get_pretty_name(population),fontsize=5,y=0.96,loc='left')
    
        axis.spines['top'].set_visible(False)
        axis.spines['right'].set_visible(False)
        axis.get_yaxis().tick_left()
        axis.get_yaxis().set_tick_params(direction='out',length=2)
        axis.set_xticks([])
        axis.set_ylabel('Fixed mutations')    
    
        axis.set_xlim([-1,22])
    
        if population_idx==5:
            axis.set_xlabel('Clones')
            
            
########################################
#
# Now do the plotting (focal first, then rest)
#
########################################

theory_times = numpy.arange(0,121)*500

gene_names, start_positions, end_positions, promoter_start_positions, promoter_end_positions, gene_sequences, strands = parse_file.parse_gene_list()

gene_name_position_map = {gene_names[i]: (start_positions[i], end_positions[i]) for i in xrange(0,len(gene_names))}

state_color_map = {parse_file.clade_hmm_states['FB']:'0.7', parse_file.clade_hmm_states['FM']:'#7a0177', parse_file.clade_hmm_states['Fm']: '#f768a1'}

for metapopulation_idx in xrange(0,2):
    for population_idx in xrange(0,6):

        population = metapopulations[metapopulation_idx][population_idx] 
        axis = axes[metapopulation_idx][population_idx]
        sys.stderr.write('Processing %s...\n' % population)

        # load mutations
        mutations, depth_tuple = parse_file.parse_annotated_timecourse(population)
    
        population_avg_depth_times, population_avg_depths, clone_avg_depth_times, clone_avg_depths = depth_tuple
    
        dummy_times,fmajors,fminors,haplotype_trajectories =  parse_file.parse_haplotype_timecourse(population)
   
        clone_mutation_Ls = []
        for clone_idx in xrange(0,len(clone_avg_depth_times)):
            clone_mutation_Ls.append( {state: 0 for state in state_color_map} )
   
        for mutation_idx in xrange(0,len(mutations)):
    
            location, gene_name, allele, var_type, test_statistic, pvalue, cutoff_idx, depth_fold_change, depth_change_pvalue, times, alts, depths, clone_times, clone_alts, clone_depths = mutations[mutation_idx] 
    
            Ls = haplotype_trajectories[mutation_idx]
    
            good_idxs, filtered_alts, filtered_depths = timecourse_utils.mask_timepoints(times, alts, depths, var_type, cutoff_idx, depth_fold_change, depth_change_pvalue)
        
            
            clone_freqs = clone_alts*1.0/(clone_depths+(clone_depths==0))
            clone_depth_fold_changes = timecourse_utils.estimate_depth_fold_changes(clone_avg_depths, clone_depths)

            for clone_idx in xrange(0,len(clone_times)):
            
                L = Ls[numpy.fabs(times-clone_times[clone_idx]).argmin()]
                
                if L not in state_color_map:
                    continue
                    
                # state is in allowed states
                # let's see whether the clone "has" it
                if clone_depths[clone_idx] > parse_file.default_min_depth:
                    if clone_depth_fold_changes[clone_idx] > -4:
                        if clone_freqs[clone_idx] > 0.5:
                            clone_mutation_Ls[clone_idx][L] += 1
        #time to plot!
        for clone_idx in xrange(0,len(clone_mutation_Ls)):
            
            x = clone_idx
            width = 0.6
            
            y0 = 0
            y1 = clone_mutation_Ls[clone_idx][parse_file.clade_hmm_states['FB']]
            color = state_color_map[parse_file.clade_hmm_states['FB']]
            
            if y1-y0 > 0.5:
                axis.fill_between([x,x+0.7],[y0,y0],[y1,y1],color=color)
            
            y0 += (y1-y0)
            y1 += clone_mutation_Ls[clone_idx][parse_file.clade_hmm_states['FM']]
            color = state_color_map[parse_file.clade_hmm_states['FM']]
            
            if y1-y0 > 0.5:
                axis.fill_between([x,x+0.7],[y0,y0],[y1,y1],color=color)
            
            y0 += (y1-y0)
            y1 += clone_mutation_Ls[clone_idx][parse_file.clade_hmm_states['Fm']]
            color = state_color_map[parse_file.clade_hmm_states['Fm']]
            
            if y1-y0 > 0.5:
                axis.fill_between([x,x+0.7],[y0,y0],[y1,y1],color=color)
        
        if metapopulation_idx==0 and population_idx==0:
            
            color = state_color_map[parse_file.clade_hmm_states['FB']]
            axis.bar([-1],[1], facecolor=color,edgecolor='none',label='$F_B$')
            
            color = state_color_map[parse_file.clade_hmm_states['FM']]
            axis.bar([-1],[1], facecolor=color,edgecolor='none',label='$F_M$')
            
            color = state_color_map[parse_file.clade_hmm_states['Fm']]
            axis.bar([-1],[1], facecolor=color,edgecolor='none',label='$F_m$')
            
            axis.legend(loc='upper left',frameon=False,fontsize=7)
            
            
sys.stderr.write("Saving final image...\t")
fig.savefig(parse_file.figure_directory+'supplemental_clone_clade_verification.pdf',bbox_inches='tight',transparent=True,dpi=300)
sys.stderr.write("Done!\n")
