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

PLOT_FMAJOR=None
PLOT_FMINOR=None
additional_titles=None
COLORED_LINEWIDTH=1

PLOT_APPEARANCE_TIME=False

# load settings
settings_filename = sys.argv[1]
settings_file = open(settings_filename,"r")
settings_string = "\n".join(settings_file.readlines())
settings_file.close()
exec settings_string    

if PLOT_FMAJOR==None:
    PLOT_FMAJOR=[False for population in populations]
if PLOT_FMINOR==None:
    PLOT_FMINOR=[False for population in populations]
if additional_titles==None:
    additional_titles = ["" for population in populations]

mpl.rcParams['font.size'] = 10.0
mpl.rcParams['lines.linewidth'] = 2.0
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'

theory_times = numpy.arange(0,121)*500

gene_names, start_positions, end_positions, promoter_start_positions, promoter_end_positions, gene_sequences, strands = parse_file.parse_gene_list()

gene_name_position_map = {gene_names[i]: (start_positions[i], end_positions[i]) for i in xrange(0,len(gene_names))}

fig_width = 6+8*(tmax-tmin)*1.0/60000
fig_height = 2*(len(populations))

fig, axes = plt.subplots(len(populations),sharex=True,sharey=True,figsize=(fig_width, fig_height))

if len(populations)<2:
    axes = [axes]

#fig = plt.figure(figsize=(fig_width,fig_height))
#outer_grid  = gridspec.GridSpec(len(populations), 1)
xticks = [5000*i for i in xrange(0,13)]
xticklabels = ['%dk' % (5*i) for i in xrange(0,13)]
xticklabels[0] = '0'
    

freq_axis = None
    
for population_idx in xrange(0,len(populations)):        

    population = populations[population_idx]

    sys.stderr.write("Processing %s...\n" % parse_file.get_pretty_name(population))

    # calculate mutation trajectories
    # load mutations
    mutations, depth_tuple = parse_file.parse_annotated_timecourse(population)
    
    population_avg_depth_times, population_avg_depths, clone_avg_depth_times, clone_avg_depths = depth_tuple
    
    dummy_times,fmajors,fminors,haplotype_trajectories = parse_file.parse_haplotype_timecourse(population)
    state_times, state_trajectories = parse_file.parse_well_mixed_state_timecourse(population)
    
    # set up figure axis
    freq_axis = axes[population_idx]
    
    if additional_titles[population_idx]=="":
        title_text = parse_file.get_pretty_name(population) 
    else:
        title_text = '%s %s' % (parse_file.get_pretty_name(population), additional_titles[population_idx])  
    freq_axis.set_title(title_text,fontsize=9,loc='right')
        
    freq_axis.set_ylabel('Allele frequency, $f(t)$',fontsize=9)
    cl = pylab.getp(freq_axis, 'ymajorticklabels')
    pylab.setp(cl, fontsize=9) 

    freq_axis.spines['top'].set_visible(False)
    freq_axis.spines['right'].set_visible(False)
    freq_axis.get_xaxis().tick_bottom()
    freq_axis.get_yaxis().tick_left()
    
    num_colored_mutations = 0
    num_total_mutations = 0

    for mutation_idx in xrange(0,len(mutations)):
     
        location, gene_name, allele, var_type, test_statistic, pvalue, cutoff_idx, depth_fold_change, depth_change_pvalue, times, alts, depths, clone_times, clone_alts, clone_depths = mutations[mutation_idx] 
        Ls = haplotype_trajectories[mutation_idx]
        state_Ls = state_trajectories[mutation_idx]
        
        good_idxs, filtered_alts, filtered_depths = timecourse_utils.mask_timepoints(times, alts, depths, var_type, cutoff_idx, depth_fold_change, depth_change_pvalue)
    
        freqs = timecourse_utils.estimate_frequencies(filtered_alts, filtered_depths)
        
        masked_times = times[good_idxs]
        masked_freqs = freqs[good_idxs]
        masked_state_Ls = state_Ls[good_idxs]
        masked_Ls = Ls[good_idxs]
        masked_depth_ratios = depths[good_idxs]/population_avg_depths[good_idxs]
            
        t = timecourse_utils.calculate_appearance_time(masked_times, masked_freqs, masked_state_Ls, masked_Ls)
        fixed_weight = timecourse_utils.calculate_fixed_weight(masked_state_Ls[-1], masked_freqs[-1])
       
        interpolation_function = timecourse_utils.create_interpolation_function(masked_times, masked_freqs, tmax=100000)
        haplotype_interpolation_function = timecourse_utils.create_interpolation_function(masked_times, masked_Ls, tmax=100000,kind='nearest')
    
        masked_clone_times, masked_clone_freqs = timecourse_utils.estimate_clone_frequencies(clone_times, clone_alts, clone_depths)
  
    
        clone_freqs = clone_alts*1.0/(clone_depths+(clone_depths==0))
        clone_depth_fold_changes = timecourse_utils.estimate_depth_fold_changes(clone_avg_depths, clone_depths)
    
        num_total_mutations +=1
      
        extra_data = (masked_clone_times, masked_clone_freqs, clone_depth_fold_changes, cutoff_idx, depth_fold_change, depth_change_pvalue, masked_freqs, masked_depth_ratios, allele)
      
        if color_condition(population_idx, location, gene_name, var_type, interpolation_function, haplotype_interpolation_function, extra_data):
        
            # One of the colored ones!
            num_colored_mutations+=1
            
            sys.stderr.write("%s %d %s %s\n" % (gene_name, location, var_type, allele)) 
            
            line, = freq_axis.plot(masked_times, masked_freqs, '-o', alpha=0.5, markersize=2, markeredgecolor='none', zorder=4, linewidth=COLORED_LINEWIDTH)
            color = pylab.getp(line,'color')
            
            if PLOT_APPEARANCE_TIME:
            
                freq_axis.plot([t],[interpolation_function(t)],'*',markersize=10,markeredgecolor='k',zorder='5',color=color)
            
            #freq_axis.plot(theory_times, interpolation_function(theory_times), '-', alpha=0.5, zorder=4, linewidth=1, color=color)
        
        else:  
            # One of the non-colored ones
            #freq_axis.plot(theory_times, interpolation_function(theory_times), '-', alpha=0.5, color='0.7', markersize=3,linewidth=1,zorder=1)
            
            freq_axis.plot(masked_times, masked_freqs, '-', color='0.7', alpha=0.5, markersize=3,label=gene_name,linewidth=1,zorder=1)
     
    sys.stderr.write("Colored=%d, Total=%d\n" % (num_colored_mutations, num_total_mutations))
    
    if PLOT_FMAJOR[population_idx]:
        freq_axis.plot(dummy_times, fmajors, 'k-', linewidth=3,zorder=10)
    if PLOT_FMINOR[population_idx]:
        freq_axis.plot(dummy_times, fminors, 'k-', linewidth=3,zorder=10)

freq_axis.set_xticks(xticks)
freq_axis.set_xticklabels(xticklabels)            
freq_axis.set_xlabel('Generation, $t$',fontsize=9)
cl = pylab.getp(freq_axis, 'xmajorticklabels')
pylab.setp(cl, fontsize=9) 

freq_axis.set_ylim([0,1.02])
freq_axis.set_xlim([tmin,tmax])   
 
sys.stderr.write("Saving final PNG image...\t")
fig.savefig(filename, bbox_inches='tight', dpi=300, transparent=True)
pylab.close(fig)
sys.stderr.write("Done!\n")

